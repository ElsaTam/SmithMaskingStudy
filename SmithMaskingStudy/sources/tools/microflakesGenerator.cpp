#include "tools/microflakesGenerator.h"

#include <random>
#include "stb_image.h"
#include "utils/math/math.h"

MicroflakesGenerator::MicroflakesGenerator(const std::string& filename)
{
	int channels;
	unsigned char* data = stbi_load(filename.c_str(), &(hf.width), &(hf.height), &channels, 1);

	hMin = std::numeric_limits<scal>::max();
	hMax = std::numeric_limits<scal>::min();

	hf.greys.resize(hf.height);
	for (int r = 0; r < hf.height; ++r) {
		hf.greys[r].resize(hf.width);
		for (int c = 0; c < hf.width; ++c) {
			int idx = c * hf.height + r;
			hf.greys[r][c] = (uint8_t)data[idx];
			if (hMin > hf.greys[r][c]) hMin = hf.greys[r][c];
			if (hMax < hf.greys[r][c]) hMax = hf.greys[r][c];
		}
	}
	stbi_image_free(data);

	size_t lastDelim = filename.find_last_of("/");
	name = filename.substr(lastDelim + 1);
	name = name.substr(0, name.find_last_of("."));
	name = name + "_microflakes";
}


template<class T>
gdt::vec3sc MicroflakesGenerator::barycentric(const T& P, const T& A, const T& B, const T& C) const
{
	gdt::vec3sc coords;
	T v0 = B - A, v1 = C - A, v2 = P - A;
	float d00 = dot(v0, v0);
	float d01 = dot(v0, v1);
	float d11 = dot(v1, v1);
	float d20 = dot(v2, v0);
	float d21 = dot(v2, v1);
	float denom = d00 * d11 - d01 * d01;
	coords.v = (d11 * d20 - d01 * d21) / denom;
	coords.w = (d00 * d21 - d01 * d20) / denom;
	coords.u = 1.0f - coords.v - coords.w;
	return coords;
}

gdt::vec3sc MicroflakesGenerator::getVertexNormal(const gdt::vec2i& pixel) const
{
	scal T = hf.scale * hf.greys[pixel.x][pixel.y - 1];
	scal B = hf.scale * hf.greys[pixel.x][pixel.y + 1];
	scal L = hf.scale * hf.greys[pixel.x - 1][pixel.y];
	scal R = hf.scale * hf.greys[pixel.x + 1][pixel.y];

	gdt::vec3sc n((L - R), (T - B), 2);
	return gdt::normalize(n);
}

gdt::vec3sc MicroflakesGenerator::getFlatNormal(const gdt::vec2sc& coords) const
{
	gdt::vec2sc pixel;
	gdt::vec2sc fractionalCoord;
	fractionalCoord.x = std::modf(coords.x, &(pixel.x));
	fractionalCoord.y = std::modf(coords.y, &(pixel.y));

	gdt::vec3sc O(pixel.x + 0.5, pixel.y + 0.5, hf.scale * hf.greys[(int)pixel.x][(int)pixel.y]);
	gdt::vec3sc T(pixel.x + 0.5, pixel.y - 0.5, hf.scale * hf.greys[(int)pixel.x][(int)pixel.y - 1]);
	gdt::vec3sc B(pixel.x + 0.5, pixel.y + 1.5, hf.scale * hf.greys[(int)pixel.x][(int)pixel.y + 1]);
	gdt::vec3sc L(pixel.x - 0.5, pixel.y + 0.5, hf.scale * hf.greys[(int)pixel.x - 1][(int)pixel.y]);
	gdt::vec3sc R(pixel.x + 1.5, pixel.y + 0.5, hf.scale * hf.greys[(int)pixel.x + 1][(int)pixel.y]);

	gdt::vec3sc n, v1, v2;

	if (fractionalCoord.x < 0.5 && fractionalCoord.y < 0.5) // bottom left
	{
		v1 = B - O;
		v2 = L - O;
	}
	else if (fractionalCoord.x < 0.5 && fractionalCoord.y >= 0.5) // top left
	{
		v1 = T - O;
		v2 = L - O;
	}
	else if (fractionalCoord.x >= 0.5 && fractionalCoord.y < 0.5) // bottom right
	{
		v1 = B - O;
		v2 = R - O;
	}
	else if (fractionalCoord.x >= 0.5 && fractionalCoord.y >= 0.5) // top right
	{
		v1 = T - O;
		v2 = R - O;
	}
	n = gdt::cross(v1, v2);
	if (n.z < 0) n = -n;
	n = gdt::normalize(n);

	return n;
}

gdt::vec3sc MicroflakesGenerator::getSmoothNormal(const gdt::vec2sc& coords) const
{
	gdt::vec2sc pixel;
	gdt::vec2sc fractionalCoord;
	fractionalCoord.x = std::modf(coords.x, &(pixel.x));
	fractionalCoord.y = std::modf(coords.y, &(pixel.y));

	gdt::vec2sc O(pixel.x + 0.5f, pixel.y + 0.5f);
	gdt::vec2sc T(pixel.x + 0.5f, pixel.y - 0.5f);
	gdt::vec2sc B(pixel.x + 0.5f, pixel.y + 1.5f);
	gdt::vec2sc L(pixel.x - 0.5f, pixel.y + 0.5f);
	gdt::vec2sc R(pixel.x + 1.5f, pixel.y + 0.5f);

	gdt::vec3sc On = getVertexNormal(gdt::vec2i((int)pixel.x, (int)pixel.y ));
	gdt::vec3sc Tn = getVertexNormal(gdt::vec2i((int)pixel.x, (int)pixel.y - 1 ));
	gdt::vec3sc Bn = getVertexNormal(gdt::vec2i((int)pixel.x, (int)pixel.y + 1 ));
	gdt::vec3sc Ln = getVertexNormal(gdt::vec2i((int)pixel.x - 1, (int)pixel.y ));
	gdt::vec3sc Rn = getVertexNormal(gdt::vec2i((int)pixel.x + 1, (int)pixel.y ));

	gdt::vec3sc n;

	if (fractionalCoord.x < 0.5 && fractionalCoord.y < 0.5) // bottom left
	{
		gdt::vec3sc baryCoords = barycentric(coords, O, B, L);
		n = baryCoords.u * On + baryCoords.v * Bn + baryCoords.w * Ln;
	}
	else if (fractionalCoord.x < 0.5 && fractionalCoord.y >= 0.5) // top left
	{
		gdt::vec3sc baryCoords = barycentric(coords, O, T, L);
		n = baryCoords.u * On + baryCoords.v * Tn + baryCoords.w * Ln;
	}
	else if (fractionalCoord.x >= 0.5 && fractionalCoord.y < 0.5) // bottom right
	{
		gdt::vec3sc baryCoords = barycentric(coords, O, B, R);
		n = baryCoords.u * On + baryCoords.v * Bn + baryCoords.w * Rn;
	}
	else if (fractionalCoord.x >= 0.5 && fractionalCoord.y >= 0.5) // top right
	{
		gdt::vec3sc baryCoords = barycentric(coords, O, T, R);
		n = baryCoords.u * On + baryCoords.v * Tn + baryCoords.w * Rn;
	}
	if (n.z < 0) n = -n;
	n = gdt::normalize(n);

	return n;
}

gdt::vec3sc MicroflakesGenerator::getPoint(const gdt::vec2sc& coords) const
{
	gdt::vec2sc pixel;
	gdt::vec2sc fractionalCoord;
	fractionalCoord.x = std::modf(coords.x, &(pixel.x));
	fractionalCoord.y = std::modf(coords.y, &(pixel.y));

	gdt::vec3sc O(pixel.x + 0.5, pixel.y + 0.5, hf.scale * hf.greys[pixel.x][pixel.y]);
	gdt::vec3sc T(pixel.x + 0.5, pixel.y - 0.5, hf.scale * hf.greys[pixel.x][pixel.y - 1]);
	gdt::vec3sc B(pixel.x + 0.5, pixel.y + 1.5, hf.scale * hf.greys[pixel.x][pixel.y + 1]);
	gdt::vec3sc L(pixel.x - 0.5, pixel.y + 0.5, hf.scale * hf.greys[pixel.x - 1][pixel.y]);
	gdt::vec3sc R(pixel.x + 1.5, pixel.y + 0.5, hf.scale * hf.greys[pixel.x + 1][pixel.y]);

	gdt::vec3sc P;

	if (fractionalCoord.x < 0.5 && fractionalCoord.y < 0.5) // bottom left
	{
		gdt::vec3sc baryCoords = barycentric(coords, { O.x, O.y }, { B.x, B.y }, { L.x, L.y });
		P = baryCoords.u * O + baryCoords.v * B + baryCoords.w * L;
	}
	else if (fractionalCoord.x < 0.5 && fractionalCoord.y >= 0.5) // top left
	{
		gdt::vec3sc baryCoords = barycentric(coords, { O.x, O.y }, { T.x, T.y }, { L.x, L.y });
		P = baryCoords.u * O + baryCoords.v * T + baryCoords.w * L;
	}
	else if (fractionalCoord.x >= 0.5 && fractionalCoord.y < 0.5) // bottom right
	{
		gdt::vec3sc baryCoords = barycentric(coords, { O.x, O.y }, { B.x, B.y }, { R.x, R.y });
		P = baryCoords.u * O + baryCoords.v * B + baryCoords.w * R;
	}
	else if (fractionalCoord.x >= 0.5 && fractionalCoord.y >= 0.5) // top right
	{
		gdt::vec3sc baryCoords = barycentric(coords, { O.x, O.y }, { T.x, T.y }, { R.x, R.y });
		P = baryCoords.u * O + baryCoords.v * T + baryCoords.w * R;
	}
	return P;
}

MicroflakesGenerator::Flake MicroflakesGenerator::createFlake(scal size, const gdt::vec3sc& P, const gdt::vec3sc& N) const
{
	// create points in the horizontal plane
	// A --- B
	// |  0  |
	// D --- C
	gdt::vec3sc A(-size * 0.5, -size * 0.5, 0);
	gdt::vec3sc B(size * 0.5,  -size * 0.5, 0);
	gdt::vec3sc C(size * 0.5,   size * 0.5, 0);
	gdt::vec3sc D(-size * 0.5,  size * 0.5, 0);

	// rotate those points with the given normal
	A = Geometry::rotateAlongNormal(A, N);
	B = Geometry::rotateAlongNormal(B, N);
	C = Geometry::rotateAlongNormal(C, N);
	D = Geometry::rotateAlongNormal(D, N);

	// center those points on P
	A = P + A;
	B = P + B;
	C = P + C;
	D = P + D;

	Flake flake{ size, N, A, B, C, D };

	return flake;
}

void MicroflakesGenerator::jiggleFlake(Flake& flake) const {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<scal> unif = std::uniform_real_distribution<scal>(hMin, hMax);

	scal hNew = unif(gen);
	scal hFlake = (flake.A.z + flake.B.z + flake.C.z + flake.D.z) / (scal)4;
	scal perturb = hNew - hFlake;
	flake.A.z += perturb;
	flake.B.z += perturb;
	flake.C.z += perturb;
	flake.D.z += perturb;
}

void MicroflakesGenerator::addFlake(const MicroflakesGenerator::Flake& flake, TriangleMesh* mesh) const
{
	// add the 4 vertices to the mesh
	// A
	mesh->vertex.push_back(flake.A);
	mesh->vertex_normal.push_back(flake.normal);
	int idxA = (int) mesh->vertex.size() - 1;
	// B
	mesh->vertex.push_back(flake.B);
	mesh->vertex_normal.push_back(flake.normal);
	int idxB = (int) mesh->vertex.size() - 1;
	// C
	mesh->vertex.push_back(flake.C);
	mesh->vertex_normal.push_back(flake.normal);
	int idxC = (int) mesh->vertex.size() - 1;
	// D
	mesh->vertex.push_back(flake.D);
	mesh->vertex_normal.push_back(flake.normal);
	int idxD = (int) mesh->vertex.size() - 1;

	// add the two triangles to the mesh
	// triangle 1
	mesh->index.push_back({ idxA, idxB, idxC });
	mesh->triangle_normal.push_back(flake.normal);
	mesh->area.push_back(flake.size * flake.size * 0.5);
	// triangle 2
	mesh->index.push_back({ idxA, idxC, idxD });
	mesh->triangle_normal.push_back(flake.normal);
	mesh->area.push_back(flake.size * flake.size * 0.5);

	// increase the total surface area
	mesh->surfaceArea += flake.size * flake.size;
}

TriangleMesh* MicroflakesGenerator::createModel(const gdt::vec2i gridSize) const
{
	// mesh initialisation
	TriangleMesh* mesh = new TriangleMesh;
	mesh->name = name;

	// find the size of a flake, we need to avoid any overlaping
	// we won't work with the border pixels
	int border = 2;
	scal stepX = (scal)(hf.width - 1 - 2 * border) / (scal)(gridSize.x - 1);
	scal stepY = (scal)(hf.height - 1 - 2 * border) / (scal)(gridSize.y - 1);
	scal flakeSize = std::min(stepX, stepY);

	// create the flakes
	for (int i = 0; i < gridSize.x; ++i)
	{
		gdt::vec2sc coords;
		coords.x = border + i * stepX;
		for (int j = 0; j < gridSize.y; ++j)
		{
			coords.y = border + j * stepY;

			gdt::vec3sc P = getPoint(coords);
			gdt::vec3sc N = getSmoothNormal(coords);
			Flake flake = createFlake(flakeSize, P, N);
			//jiggleFlake(flake);
			addFlake(flake, mesh);
		}
	}

	mesh->surfaceArea = flakeSize * gridSize.x * gridSize.y;

	// compute the mesonormal
	mesh->meso_normal = fittingPlaneNormal(mesh->vertex);

	// extend the bounding box
	for (auto vtx : mesh->vertex)
		mesh->bounds.extend(vtx);

	return mesh;
}