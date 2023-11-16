#pragma once

#include <string>
#include <vector>
#include "shapes/TriangleMesh.h"

class MicroflakesGenerator
{
private:

	struct Heightfield {
		int width = 0;
		int height = 0;
		std::vector<std::vector<uint8_t>> greys;
		scal scale = 0.07808f;
	} hf;

	scal hMin;
	scal hMax;

	std::string name;

	struct Flake {
		scal size;
		gdt::vec3sc normal;
		gdt::vec3sc A;
		gdt::vec3sc B;
		gdt::vec3sc C;
		gdt::vec3sc D;
	};

	/**
	 * @brief    Compute barycentric coordinates (u, v, w) for point P with respect to triangle (A, B, C)
	 * @tparam T Vector in a specific dimension (vec2sc, vec3sc, ...)
	 * @param P  The point we want the coordinates of
	 * @param A  First point of the triangle
	 * @param B  Second point of the triangle
	 * @param C  Third point of the triangle
	 * @return   Barycentric coordinates (u, v, w)
	*/
	template<class T>
	gdt::vec3sc barycentric(const T& P, const T& A, const T& B, const T& C) const;

	/**
	 * @brief Compute the normal for a vertex (at the center of a pixel)
	 * @param pixel 
	 * @return 
	*/
	gdt::vec3sc getVertexNormal(const gdt::vec2i& pixel) const;

	/**
	 * @brief Compute the flat normal at the sampled coordinates. Computed with 1 face normal.
	 * @param coords 
	 * @return The normal of the triangle that contains (coords, h(coords))
	*/
	gdt::vec3sc getFlatNormal(const gdt::vec2sc& coords) const;

	/**
	 * @brief Compute the smooth normal at a point. Computed with 3 vertex normals.
	*/
	gdt::vec3sc getSmoothNormal(const gdt::vec2sc& coords) const;

	/**
	 * @brief Get the 3D point at the sampled coordinates.
	*/
	gdt::vec3sc getPoint(const gdt::vec2sc& coords) const;

	/**
	 * @brief Create a flake at point P, with normal N
	*/
	Flake createFlake(scal size, const gdt::vec3sc& P, const gdt::vec3sc& N) const;

	void jiggleFlake(Flake& flake) const;

	/**
	 * @brief Add a flake (two triangles and their normal) to a mesh
	*/
	void addFlake(const Flake& flake, TriangleMesh* mesh) const;

public:
	MicroflakesGenerator(const std::string& filename);

	TriangleMesh* createModel(const gdt::vec2i gridSize) const;
};