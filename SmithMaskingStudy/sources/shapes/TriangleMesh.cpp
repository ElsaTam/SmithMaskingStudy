// ======================================================================== //
// Copyright 2018-2019 Ingo Wald                                            //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include <random>
#include <numeric>
#include <unordered_set>
#include <Eigen/Dense>
#include <filesystem>

#include "shapes/TriangleMesh.h"
#include "utils/console.h"
#include "utils/params.h"
#include "utils/kMeans.h"
#include "tools/logger.h"
#include "tools/microflakesGenerator.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

#include "igl/readOBJ.h"
#include "igl/per_face_normals.h"
#include "igl/per_vertex_normals.h"
#include "igl/doublearea.h"

using namespace gdt;


namespace std {
    inline bool operator<(const tinyobj::index_t &a,
                        const tinyobj::index_t &b)
    {
        if (a.vertex_index < b.vertex_index) return true;
        if (a.vertex_index > b.vertex_index) return false;
    
        if (a.normal_index < b.normal_index) return true;
        if (a.normal_index > b.normal_index) return false;
    
        if (a.texcoord_index < b.texcoord_index) return true;
        if (a.texcoord_index > b.texcoord_index) return false;
    
        return false;
    }
}


/*! find vertex with given position, normal, texcoord, and return
    its vertex ID, or, if it doesn't exit, add it to the mesh, and
    its just-created index */
int addVertex(TriangleMesh *mesh,
            tinyobj::attrib_t &attributes,
            const tinyobj::index_t &idx,
            const tinyobj::index_t &otherIdx1,
            const tinyobj::index_t &otherIdx2,
            std::map<int,int> &knownVertices)
{
    const vec3sc *vertex_array   = (const vec3sc*)attributes.vertices.data();
    const vec3sc *normal_array   = (const vec3sc*)attributes.normals.data();

    if (knownVertices.find(idx.vertex_index) != knownVertices.end()) { // vertex already in the set
        int idxInMesh = knownVertices[idx.vertex_index];
        if (idx.normal_index >= 0) {
            vec3sc n = normalize(normal_array[idx.normal_index]);
            if (n.z < 0) n = -n;
            mesh->vertex_normal[idxInMesh] += n;
        }
        else {
            const vec3sc& A = vertex_array[idx.vertex_index];
            const vec3sc& B = vertex_array[otherIdx1.vertex_index];
            const vec3sc& C = vertex_array[otherIdx2.vertex_index];
            vec3sc sideA = B - A, sideB = C - A;
            vec3sc n = normalize(cross(sideA, sideB));
            if (n.z < 0) n = -n;
            mesh->vertex_normal[idxInMesh] += n;
        }
        return idxInMesh;
    }
    
    int newID = (int) mesh->vertex.size();
    knownVertices[idx.vertex_index] = newID;

    mesh->vertex.push_back(vertex_array[idx.vertex_index]);

    if (idx.normal_index >= 0) {
        vec3sc n = normalize(normal_array[idx.normal_index]);
        if (n.z < 0) n = -n;
        mesh->vertex_normal.push_back(n);
    }
    else {
        const vec3sc& A = vertex_array[idx.vertex_index];
        const vec3sc& B = vertex_array[otherIdx1.vertex_index];
        const vec3sc& C = vertex_array[otherIdx2.vertex_index];
        vec3sc sideA = B - A, sideB = C - A;
        vec3sc n = normalize(cross(sideA, sideB));
        if (n.z < 0) n = -n;
        mesh->vertex_normal.push_back(n);
    }

    // just for sanity's sake:
    if (mesh->vertex_normal.size() > 0)
        mesh->vertex_normal.resize(mesh->vertex.size());
    
    return newID;
}

scal computeArea(TriangleMesh* mesh, vec3i idx)
{
    const vec3sc& A = mesh->vertex[idx[0]];
    const vec3sc& B = mesh->vertex[idx[1]];
    const vec3sc& C = mesh->vertex[idx[2]];
    vec3sc sideA = B - A, sideB = C - A;
    return 0.5 * length(cross(sideA, sideB));
}

vec3sc triangleNormal(TriangleMesh* mesh, vec3i idx)
{
    const vec3sc& A = mesh->vertex[idx[0]];
    const vec3sc& B = mesh->vertex[idx[1]];
    const vec3sc& C = mesh->vertex[idx[2]];
    vec3sc sideA = B - A, sideB = C - A;
    vec3sc N = normalize(cross(sideA, sideB));
    if (dot(N, mesh->meso_normal) < 0) N = -N;
    return N;
}

vec3sc fittingPlaneNormal(const std::vector<vec3sc>& vertices)
{
    //return { 0, 0, 1 };
    // copy coordinates to  matrix in Eigen format
    size_t N = vertices.size();
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > coord(3, N);
    for (size_t i = 0; i < N; ++i)
        coord.col(i) = Eigen::Vector3d(vertices[i].x, vertices[i].y, vertices[i].z);

    // calculate centroid
    Eigen::Vector3d centroid(coord.row(0).mean(), coord.row(1).mean(), coord.row(2).mean());

    // subtract centroid
    coord.row(0).array() -= centroid(0);
    coord.row(1).array() -= centroid(1);
    coord.row(2).array() -= centroid(2);

    // we only need the left-singular matrix here
    //  http://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
    auto svd = coord.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::Vector3d n = svd.matrixU().rightCols<1>();
    if (n(2) < 0) n = -n;

    vec3sc normal((scal)n(0), (scal)n(1), (scal)n(2));
    return normal;
}

vec3sc meanNormal(const std::vector<vec3sc>& normals, const std::vector<scal>& areas) {
    vec3sc N(0);
    scal A = 0;
    for (int i = 0; i < normals.size(); ++i) {
        N += normals[i] * areas[i];
    }
    N /= normals.size();
    N = normalize(N);
    return N;
}

vec3sc TriangleMesh::computeMesoNormal(int maxSamples) const {
    std::vector<vec3sc> pointsCloud;
    for (const gdt::vec3sc& v : vertex) {
        if (bounds.closest_distance(v).x < Parameters::get()->currentParams()->sideEffectParams.borderPercentage * bounds.span().x / 2.f
            || bounds.closest_distance(v).y < Parameters::get()->currentParams()->sideEffectParams.borderPercentage * bounds.span().y / 2.f)
        {
            continue;
        }
        pointsCloud.push_back(v);
    }
    if (maxSamples > 0) {
        std::default_random_engine rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<scal> rng(0.f, 1.f);
        for (int faceID = 0; faceID < this->index.size(); ++faceID) {
            const vec3i& idx = this->index[faceID];
            const scal A = this->area[faceID];
            /* Uniform sampling from unit-square */
            const vec3sc& v1 = this->vertex[idx[0]];
            const vec3sc& v2 = this->vertex[idx[1]];
            const vec3sc& v3 = this->vertex[idx[2]];
            int nSamples = (maxSamples * A) / this->maxArea;
            for (int i = 0; i < nSamples; ++i) {
                scal s1 = rng(gen);
                scal s2 = rng(gen);
                /* Mapping to unit-triangle */
                scal t1 = 0.5f * s1;
                scal t2 = 0.5f * s2;
                scal offset = t2 - t1;
                if (offset > 0) { t2 += offset; }
                else { t1 -= offset; }
                /* Mapping to arbitrary triangle */
                vec3sc point(t1 * v1 + t2 * v2 + (1.f - t1 - t2) * v3);
                pointsCloud.push_back(point);
            }
        }
    }
    return fittingPlaneNormal(pointsCloud);
}

TriangleMesh* loadOBJ(const std::string &objFile)
{
    const auto start = std::chrono::high_resolution_clock::now();

    const bool flip = false;

    TriangleMesh *_mesh = new TriangleMesh;

    std::filesystem::path p = objFile;
    _mesh->name = p.stem().string();
    if (flip) {
        _mesh->name += "_flip";
    }

    tinyobj::attrib_t attributes;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string err = "";

    Console::print(OutLevel::TRACE, Console::timeStamp.str() + "tinyobj::LoadObj...");
    bool readOK
        = tinyobj::LoadObj(&attributes,
                           &shapes,
                           &materials,
                           &err,
                           &err,
						   objFile.c_str(),
                           NULL,
                           true);
    if (!readOK) {
        throw std::runtime_error("Could not read OBJ model from " + objFile + " : " + err);
    }

    // flip surface
    if (flip) {
        Console::print(OutLevel::TRACE, Console::timeStamp.str() + "Flip surface...");
        for (int vID = 0; vID < attributes.vertices.size() / 3; vID++) {
            attributes.vertices[(size_t)(3 * vID + 2)] *= -1; // flip Z
        }
    }
    
    Console::print(OutLevel::TRACE, Console::timePad + "[tinyobj] shape.size()  = " + std::to_string(shapes.size()));
    Console::print(OutLevel::TRACE, Console::timePad + "[tinyobj] attributes.vertices.size() = " + std::to_string(attributes.vertices.size()) + " (" + std::to_string(attributes.vertices.size() / 3) + " * 3)");
    Console::print(OutLevel::TRACE, Console::timePad + "[tinyobj] attributes.normals.size()  = " + std::to_string(attributes.normals.size()) + " (" + std::to_string(attributes.normals.size() / 3) + " * 3)");
    Console::print(OutLevel::TRACE, Console::timePad + "[tinyobj] shape.mesh.indices.size()  = " + std::to_string(shapes[0].mesh.indices.size()) + " (" + std::to_string(shapes[0].mesh.indices.size() / 3) + " * 3)");

    tinyobj::shape_t& shape = shapes[0]; // only one shape
      
    std::map<int,int> knownVertices;

    _mesh->surfaceArea = 0;
    _mesh->macroArea = 0;
    _mesh->maxArea = 0;

    // create mesh
    Console::print(OutLevel::TRACE, Console::timeStamp.str() + "Add vertices...");
    for (int faceID = 0; faceID < shape.mesh.indices.size() / 3; faceID++) {
        tinyobj::index_t idx0 = shape.mesh.indices[(size_t)(3*faceID+0)];
        tinyobj::index_t idx1 = shape.mesh.indices[(size_t)(3*faceID+1)];
        tinyobj::index_t idx2 = shape.mesh.indices[(size_t)(3*faceID+2)];
          
        vec3i idx(addVertex(_mesh, attributes, idx0, idx1, idx2, knownVertices),
                  addVertex(_mesh, attributes, idx1, idx2, idx0, knownVertices),
                  addVertex(_mesh, attributes, idx2, idx0, idx1, knownVertices));

        scal A = computeArea(_mesh, idx);
        _mesh->surfaceArea += A;
        _mesh->area.push_back(A);
        _mesh->index.push_back(idx);

        if (_mesh->maxArea < A) _mesh->maxArea = A;
    }

    // normalize each vertex normal
    Console::print(OutLevel::TRACE, Console::timeStamp.str() + "Normalize vertex normals...");
    for (int normalID = 0; normalID < _mesh->vertex_normal.size(); ++normalID) {
        _mesh->vertex_normal[normalID] = normalize(_mesh->vertex_normal[normalID]);
    }

    // find the best fitting plane (and its mesonormal)
    Console::print(OutLevel::TRACE, Console::timeStamp.str() + "Mesonormal...");
    _mesh->meso_normal = _mesh->computeMesoNormal(0);


    // compute normals for each triangle
    // must be done after the mesonormal, because it is used to choose the sign of the micronormals.
    Console::print(OutLevel::TRACE, Console::timeStamp.str() + "Triangle normals...");
    for (int faceID = 0; faceID < _mesh->index.size(); faceID++) {
        _mesh->triangle_normal.push_back(triangleNormal(_mesh, _mesh->index[faceID]));
    }


    // assign the mesh
    if (_mesh->vertex.empty()) {
        delete _mesh;
        if (Parameters::get()->currentParams()->outLevel >= OutLevel::ERR)
            Console::print(OutLevel::ERR,"Done loading obj file - no vertex found.");
        return nullptr;
    }

    // extend the bounding box
    Console::print(OutLevel::TRACE, Console::timeStamp.str() + "BBox...");
    for (auto vtx : _mesh->vertex)
        _mesh->bounds.extend(vtx);
    _mesh->macroArea = _mesh->bounds.span().x * _mesh->bounds.span().y;

    //model->createPatches();

    LOG_TIME_OBJ(_mesh->name + ".obj", start, std::chrono::high_resolution_clock::now());

    return _mesh;
}

TriangleMesh* createSubMesh(const TriangleMesh* _mesh, std::set<int> faceIDs)
{
    TriangleMesh* newMesh = new TriangleMesh;

    std::map<int, int> permutation; // first  --> index in _mesh
                                    // second --> index in newMesh

    newMesh->maxArea = 0;
    newMesh->surfaceArea = 0;
    newMesh->macroArea = 0;

    for (int faceID : faceIDs) {
        const vec3i& face = _mesh->index[faceID];
        vec3i newFace;
        // add the 3 vertices and their normal
        for (int i = 0; i < 3; ++i) {
            bool alreadyUsed = false;
            auto itIndex = permutation.find(face[i]);
            if (itIndex != permutation.end()) {
                newFace[i] = permutation[face[i]];
            }
            else {
                newMesh->vertex.push_back(_mesh->vertex[face[i]]);
                newMesh->vertex_normal.push_back(_mesh->vertex_normal[face[i]]);
                newFace[i] = (int)newMesh->vertex.size() - 1; // update the face with correct index
                permutation[face[i]] = newFace[i];
            }
        }
        newMesh->index.push_back(newFace); // add the face indices
        newMesh->triangle_normal.push_back(_mesh->triangle_normal[faceID]); // add the face normal
        newMesh->area.push_back(_mesh->area[faceID]); // add the face area
        newMesh->surfaceArea += _mesh->area[faceID];
        if (_mesh->area[faceID] > newMesh->maxArea)
            newMesh->maxArea = _mesh->area[faceID];
    }
    newMesh->meso_normal = fittingPlaneNormal(newMesh->vertex);

    newMesh->name = _mesh->name;

    // extend the bounding box
    for (auto vtx : newMesh->vertex)
        newMesh->bounds.extend(vtx);
    newMesh->macroArea = newMesh->bounds.span().x * newMesh->bounds.span().y;

    Console::print(OutLevel::INFO, Console::timeStamp.str() + "Old mesh: ");
    Console::print(OutLevel::INFO, Console::timeStamp.str() + "Vertex: " + std::to_string(_mesh->vertex.size()));
    Console::print(OutLevel::INFO, Console::timeStamp.str() + "Faces: " + std::to_string(_mesh->index.size()));
    std::stringstream ss_bounds;
    ss_bounds << _mesh->bounds;
    Console::print(OutLevel::INFO, Console::timeStamp.str() + "Bounds: " + ss_bounds.str());
    Console::print(OutLevel::INFO, Console::timeStamp.str() + "New mesh: ");
    Console::print(OutLevel::INFO, Console::timeStamp.str() + "Vertex: " + std::to_string(newMesh->vertex.size()));
    Console::print(OutLevel::INFO, Console::timeStamp.str() + "Faces: " + std::to_string(newMesh->index.size()));
    ss_bounds.str("");
    ss_bounds << newMesh->bounds;
    Console::print(OutLevel::INFO, Console::timeStamp.str() + "Bounds: " + ss_bounds.str());

    return newMesh;
}


void TriangleMesh::scale(gdt::vec3sc factor) {
    vec3sc center = bounds.center();
    bounds = box3sc();
    for (size_t i = 0; i < vertex.size(); ++i) {
        vec3sc v = vertex[i];
        vec3sc vec = v - center;
        vec3sc vec_scaled = vec * factor;
        vec3sc v_new = center + vec_scaled;
        vertex[i] = center + (vertex[i] - center) * factor;
        bounds.extend(vertex[i]);
    }
    macroArea = bounds.span().x * bounds.span().y;
    surfaceArea = 0;
    area.clear();
    triangle_normal.clear();
    for (int faceID = 0; faceID < index.size(); ++faceID) {
        triangle_normal.push_back(triangleNormal(this, index[faceID]));
        area.push_back(computeArea(this, index[faceID]));
        surfaceArea += area[area.size() - 1];
    }
    meso_normal = computeMesoNormal(0);
}

void TriangleMesh::translate(vec3sc vector) {
    bounds = box3sc();
    for (size_t i = 0; i < vertex.size(); ++i) {
        vertex[i] = vertex[i] + vector;
        bounds.extend(vertex[i]);
    }
}

void TriangleMesh::rotate(vec3sc newNormal) {
    if (meso_normal == newNormal) return;

    gdt::vec3sc v = cross(meso_normal, newNormal);
    scal s = length(v);
    scal c = dot(meso_normal, newNormal);

    Eigen::Matrix<scal, 3, 3> V;
    V << 0, -v.z, v.y,
        v.z, 0, -v.x,
        -v.y, v.x, 0;

    Eigen::Matrix<scal, 3, 3> I = Eigen::Matrix<scal, 3, 3>::Identity();

    Eigen::Matrix<scal, 3, 3> R;
    R = I + V + V * V * (1. / (1. + c));

    vec3sc center = bounds.center();
    bounds = box3sc();
    for (size_t i = 0; i < vertex.size(); ++i) {
        Eigen::Matrix<scal, 3, 1> originalVertex;
        originalVertex << vertex[i].x - center.x, vertex[i].y - center.y, vertex[i].z - center.z;
        Eigen::Matrix<scal, 3, 1> orientedVertex;
        orientedVertex = R * originalVertex;
        vertex[i] = center + vec3sc(orientedVertex(0), orientedVertex(1), orientedVertex(2));

        Eigen::Matrix<scal, 3, 1> originalNormal;
        originalNormal << vertex_normal[i].x, vertex_normal[i].y, vertex_normal[i].z;
        Eigen::Matrix<scal, 3, 1> orientedNormal;
        orientedNormal = R * originalNormal;
        vertex_normal[i] = vec3sc(orientedNormal(0), orientedNormal(1), orientedNormal(2));

        bounds.extend(vertex[i]);
    }
    macroArea = bounds.span().x * bounds.span().y;
    for (size_t i = 0; i < triangle_normal.size(); ++i) {
        Eigen::Matrix<scal, 3, 1> originalNormal;
        originalNormal << triangle_normal[i].x, triangle_normal[i].y, triangle_normal[i].z;
        Eigen::Matrix<scal, 3, 1> orientedNormal;
        orientedNormal = R * originalNormal;
        triangle_normal[i] = vec3sc(orientedNormal(0), orientedNormal(1), orientedNormal(2));
    }

    meso_normal = this->computeMesoNormal(20);
}


int TriangleMesh::subdivisions() const
{
    return log2(sqrt(index.size() * 0.5));
}

void TriangleMesh::createPatches()
{
    int N = 10; // 10 x 10 patches;

    // create bounding boxes
    scal dWidth = bounds.span().x / (scal)10;
    scal dHeight = bounds.span().y / (scal)10;
    std::vector<box3sc> boxes;
    scal zmax = bounds.upper.z;
    scal zmin = bounds.lower.z;
    for (int col = 0; col < N; ++col)
    {
        for (int row = 0; row < N; ++row)
        {
            box3sc box;
            vec2sc top_left(    bounds.lower.x + dWidth * row,         bounds.lower.y + dHeight * (col + 1));
            vec2sc top_right(   bounds.lower.x + dWidth * (row + 1),   bounds.lower.y + dHeight * (col + 1));
            vec2sc bot_left(    bounds.lower.x + dWidth * row,         bounds.lower.y + dHeight * col);
            vec2sc bot_right(   bounds.lower.x + dWidth * (row + 1),   bounds.lower.y + dHeight * col);
            box.extend(vec3sc(top_left.x, top_left.y, zmax));
            box.extend(vec3sc(top_left.x, top_left.y, zmin));
            box.extend(vec3sc(top_right.x, top_right.y, zmax));
            box.extend(vec3sc(top_right.x, top_right.y, zmin));
            box.extend(vec3sc(bot_left.x, bot_left.y, zmax));
            box.extend(vec3sc(bot_left.x, bot_left.y, zmin));
            box.extend(vec3sc(bot_right.x, bot_right.y, zmax));
            box.extend(vec3sc(bot_right.x, bot_right.y, zmin));
            boxes.push_back(box);
        }
    }

    // gather vertices of each patch
    std::vector< std::vector<vec3sc> > verticesBundles(N * N);
    for (vec3sc& v : vertex)
    {
        for (int b = 0; b < boxes.size(); ++b)
        {
            if (boxes[b].contains(v)) {
                verticesBundles[b].push_back(v);
            }
        }
    }

    // update model
    std::vector< std::vector<vec3sc> > normals(N * N);
    for (int p = 0; p < N * N; ++p)
    {
        vec3sc wn = fittingPlaneNormal(verticesBundles[p]);
        patches.push_back(std::make_pair(boxes[p], wn));
    }
}

// Check if two faces are neighbor
// Parameters are the three indices of the face's vertices
bool areNeighbors(const TriangleMesh* mesh, vec3i faceA, vec3i faceB)
{
    vec3sc vA0 = mesh->vertex[faceA[0]];
    vec3sc vA1 = mesh->vertex[faceA[1]];
    vec3sc vA2 = mesh->vertex[faceA[2]];
    vec3sc vB0 = mesh->vertex[faceB[0]];
    vec3sc vB1 = mesh->vertex[faceB[1]];
    vec3sc vB2 = mesh->vertex[faceB[2]];

    int sharedV = 0;
    if (vA0 == vB0 || vA1 == vB0 || vA2 == vB0) {
        ++sharedV;
    }
    if (vA0 == vB1 || vA1 == vB1 || vA2 == vB1) {
        ++sharedV;
    }
    if (vA0 == vB2 || vA1 == vB2 || vA2 == vB2) {
        ++sharedV;
    }

    return sharedV >= 2;
}

bool areSameFlatSet(const TriangleMesh* mesh, int faceA, int faceB)
{
    const vec3sc& nA = mesh->triangle_normal[faceA];
    const vec3sc& nB = mesh->triangle_normal[faceB];
    scal d = abs(dot(nA, nB));
    return d > cos(20 * m_pi / 180);
}

bool belongsToFlatSet(const TriangleMesh* mesh, vec3sc N, int faceID)
{
    scal d = abs(dot(mesh->triangle_normal[faceID], N));
    return d > cos(80 * m_pi / 180);
}


// find the three (at most) neighbor faces
// Return the faces indices in mesh->index
std::vector<int> neighborsOfFace(const TriangleMesh* mesh, int faceID)
{
    std::vector<int> neighborFaces;

    for (int i = 0; i < mesh->index.size(); ++i) {
        if (i != faceID) {
            if (areNeighbors(mesh, mesh->index[faceID], mesh->index[i])) {
                neighborFaces.push_back(i);
            }
        }
    }

    return neighborFaces;
}

void TriangleMesh::flatSetContainingFace(int faceID, vec3sc N, std::queue<int>& faceIDQueue) const
{
    // trouver les voisins du triangle
    std::vector<int> neighbors = neighborsOfFace(this, faceID);
    for (std::vector<int>::iterator it = neighbors.begin(); it != neighbors.end(); ){
        if (!areSameFlatSet(this, faceID, *it)) {
            it = neighbors.erase(it);
        }
        else {
            faceIDQueue.push(*it);
            ++it;
        }
    }
}

template<typename T> bool contains(std::vector<T>& v, T value) {
    return std::find(v.begin(), v.end(), value) != v.end();
}
template<typename T> bool contains(std::set<T>& v, T value) {
    return v.find(value) != v.end();
}
template<typename T> void removeValue(std::vector<T>& v, T value) {
    v.erase(std::remove(v.begin(), v.end(), value), v.end());
}

std::vector<std::vector<int>> TriangleMesh::flatSets(vec3sc N) const
{
    std::vector<std::vector<int>> sets;
    std::vector<int> faces(index.size());
    std::iota(std::begin(faces), std::end(faces), 0);

    // pour chaque face
    for (std::vector<int>::iterator itFace = faces.begin(); itFace != faces.end(); ) {
        int baseFaceID = *itFace;
        // on verifie si la face est deja dans un set
        bool faceFound = false;
        for (std::vector<int>& set : sets) {
            if (contains(set, baseFaceID)) {
                faceFound = true;
                break;
            }
        }
        if (faceFound)
            continue;

        std::vector<int> flatSet;
        std::queue<int> queueSet;
        queueSet.push(baseFaceID);
        while (!queueSet.empty()) {
            int currentFaceID = queueSet.front();
            // La face a traiter est la premiere de la file.
            if (contains(flatSet, currentFaceID)) {
                // Si elle est deja dans le set, on ne fait rien
            }
            // sinon, on la traite
            else {
                // on l'ajoute au set
                flatSet.push_back(currentFaceID);
                // on cherche ses voisins
                std::vector<int> neighbors = neighborsOfFace(this, currentFaceID);
                for (int n : neighbors) {
                    // si le voisin est deja dans le set, on passe
                    if (contains(flatSet, n)) {
                        continue;
                    }
                    // on ajoute ceux qui ne sont pas plats
                    if (areSameFlatSet(this, currentFaceID, n)) {
                        queueSet.push(n);
                        // et on les enleve de la liste des faces
                        removeValue(faces, n);
                    }
                }
            }
            // on enleve la face a traiter de la file
            queueSet.pop();
        }

        // si un set a ete consitue, on l'ajoute au resultat final
        if (flatSet.size() > 10) {
            sets.push_back(flatSet);
        }

        // on passe à la face suivante
        itFace++;
    }

    return sets;
}


std::vector<std::set<int>> TriangleMesh::heightSeparation(int K, std::vector<scal>& centroids) const
{
    std::vector<std::set<int>> facesCluster(K);

    const int dim = 1;
    std::vector< vec_t<scal, dim> > heights;
    interval<int> extremaPointIDs;
    interval<scal> extremaHeights(bounds.upper.z - 1, bounds.lower.z);
    for (vec3sc p : vertex) {
        heights.push_back({ p.z });
        if (p.z < extremaHeights.lower) {
            extremaHeights.lower = p.z;
            extremaPointIDs.lower = (int)heights.size() - 1;
        }
        else if (p.z > extremaHeights.upper) {
            extremaHeights.upper = p.z;
            extremaPointIDs.upper = (int)heights.size() - 1;
        }
    }

    int iterations = 100;
    KMeans<dim> kmeans(K, iterations);
    if (K == 2) {
        std::vector<int> centersID { extremaPointIDs.lower, extremaPointIDs.upper };
        kmeans.run(heights, centersID);
    }
    else {
        kmeans.run(heights);
    }

    std::vector<std::set<int>> clusters(K);
    for (int k = 0; k < K; ++k) {
        clusters[k] = kmeans.getPointsIdOfCluster(k);
        centroids[k] = kmeans.getCentroidOfCluster(k)[0];
        Console::print(OutLevel::TRACE, Console::timePad + "Cluster " + std::to_string(k + 1) + " : " + std::to_string(clusters[k].size()) + " sommets / Centroid : " + std::to_string(centroids[k]));
    }

    int faceID;
    const size_t nFace = index.size();

#pragma omp parallel num_threads(8)
    {
        std::vector<std::set<int>> facesCluster_private(K);
        size_t local_count = 0;

        #pragma omp for nowait
        for (faceID = 0; faceID < nFace; ++faceID) {
            const vec3i& face = index[faceID];
            const scal h0 = vertex[face[0]].z, h1 = vertex[face[1]].z, h2 = vertex[face[2]].z;

            int cluster = -1;
            // le cluster auquel appartient le premier sommet
            for (int k = 0; k < K; ++k) {
                if (contains(clusters[k], face[0])) {
                    cluster = k;
                    break;
                }
            }
            // on verifie que les deux autres sommets appartiennent au meme cluster
            if (contains(clusters[cluster], face[1]) && contains(clusters[cluster], face[2])) {
                facesCluster_private[cluster].insert(faceID);
            }
        }

        #pragma omp critical
        for (int k = 0; k < K; ++k) {
            facesCluster[k].insert(facesCluster_private[k].begin(), facesCluster_private[k].end());
        }
    }

    return facesCluster;
}



vec3sc TriangleMesh::neighborhoodMesoNormal(int faceID) const
{
    const vec3i& _index = index[faceID];
    const vec3sc& A = vertex[_index[0]];
    const vec3sc& B = vertex[_index[1]];
    const vec3sc& C = vertex[_index[2]];
    const vec3sc& barycenter = (A + B + C) / (scal)3;
    for (int p = 0; p < patches.size(); ++p)
    {
        if (patches[p].first.contains(barycenter))
            return patches[p].second;
    }
    return { 0, 0, 1 };
}


vec3sc TriangleMesh::neighborhoodMesoNormal(int faceID, float d) const
{
    const vec3sc barycenter = (vertex[index[faceID][0]] + vertex[index[faceID][1]] + vertex[index[faceID][2]]) / (scal)3;
    std::unordered_set<vec3sc> patchVertices_set;
    for (int i = 0; i < index.size(); ++i)
    {
        const vec3sc A = vertex[index[i][0]];
        const vec3sc B = vertex[index[i][1]];
        const vec3sc C = vertex[index[i][2]];
        const vec3sc p = (A + B + C) / (scal)3;
        if (length(barycenter - p) <= d) {
            patchVertices_set.insert(A);
            patchVertices_set.insert(B);
            patchVertices_set.insert(C);
        }
    }
    std::vector<vec3sc> patchVertices;
    patchVertices.reserve(patchVertices_set.size());
    for (auto it = patchVertices_set.begin(); it != patchVertices_set.end(); ) {
        patchVertices.push_back(std::move(patchVertices_set.extract(it++).value()));
    }
    return fittingPlaneNormal(patchVertices);
}


vec3sc sampleTriangle(TriangleMesh* mesh, const vec3i& idx, const vec2sc& seed, vec3sc& n)
{
    const vec3sc& p0 = mesh->vertex[idx[0]];
    const vec3sc& p1 = mesh->vertex[idx[1]];
    const vec3sc& p2 = mesh->vertex[idx[2]];

    scal a = sqrt(std::max(0.0, 1.0 - seed.x));
    vec2sc bary(1 - a, a * seed.y);
    vec3sc sideA = p1 - p0, sideB = p2 - p0;
    vec3sc p = p0 + (sideA * bary.x) + (sideB * bary.y);

    if (mesh->vertex_normal.size() > 0) {
        const vec3sc& n0 = mesh->vertex_normal[idx[0]];
        const vec3sc& n1 = mesh->vertex_normal[idx[1]];
        const vec3sc& n2 = mesh->vertex_normal[idx[2]];

        n = normalize(vec3sc(
            n0 * (1.0f - bary.x - bary.y) +
            n1 * bary.x + n2 * bary.y
        ));
    }
    else {
        n = normalize(cross(sideA, sideB));
    }

    return p;
}


#include "tools/statistics.h"
TriangleMesh* createMesh(const std::string& objPath)
{
    TriangleMesh* mesh;

    // Else, building from .obj file
    Console::print(OutLevel::TRACE, Console::timeStamp.str() + "Building triangle mesh...");

    mesh = loadOBJ(objPath);


    /*
    std::vector<std::pair<std::string, scal>> scales{
        std::make_pair<std::string, scal>("Rock017_4K_Displacement", 0.2),
        std::make_pair<std::string, scal>("Rock023_4K_Displacement", 0.2),
        std::make_pair<std::string, scal>("Bark010_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Bricks078_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Bricks084_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Ground048_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Ground054_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("PavingStones055_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Rock018_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Rock022_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Rock024_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Rock025_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Rock027_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Rock040_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Rock042L_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Rock043S_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Rock045_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Rock050_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Wood086_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Ground048_4K_Displacement", 0.3),
        std::make_pair<std::string, scal>("Bark007_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("Bark009_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("Bark012_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("Bricks021_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("Ground025_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("Ground029_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("Ground052_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("Ground056_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("Ground059_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("PavingStones115A_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("PavingStones121_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("Planks028_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("Rocks023_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("Tiles039_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("TreeEnd004_4K_Displacement", 0.4),
        std::make_pair<std::string, scal>("Bark005_4K_Displacement", 0.5),
        std::make_pair<std::string, scal>("Ground038_4K_Displacement", 0.5),
        std::make_pair<std::string, scal>("PavingStones053_4K_Displacement", 0.5),
        std::make_pair<std::string, scal>("PavingStones070_4K_Displacement", 0.5),
        std::make_pair<std::string, scal>("Rocks006_4K_Displacement", 0.5),
        std::make_pair<std::string, scal>("Tiles026_4K_Displacement", 0.5),
        std::make_pair<std::string, scal>("Tiles115_4K_Displacement", 0.5),
        std::make_pair<std::string, scal>("Tiles116_4K_Displacement", 0.5),
        std::make_pair<std::string, scal>("Rock043S_4K_Displacement", 0.5),
        std::make_pair<std::string, scal>("Bricks022_4K_Displacement", 0.6),
        std::make_pair<std::string, scal>("Bricks051_4K_Displacement", 0.6),
        std::make_pair<std::string, scal>("Concrete013_4K_Displacement", 0.6),
        std::make_pair<std::string, scal>("Ground039_4K_Displacement", 0.6),
        std::make_pair<std::string, scal>("Ground055S_4K_Displacement", 0.6),
        std::make_pair<std::string, scal>("PavingStones021_4K_Displacement", 0.6),
        std::make_pair<std::string, scal>("PavingStones027_4K_Displacement", 0.6),
        std::make_pair<std::string, scal>("PavingStones044_4K_Displacement", 0.6),
        std::make_pair<std::string, scal>("PavingStones069_4K_Displacement", 0.6),
        std::make_pair<std::string, scal>("PavingStones122_4K_Displacement", 0.6),
        std::make_pair<std::string, scal>("Planks021_4K_Displacement", 0.6),
        std::make_pair<std::string, scal>("Rock046L_4K_Displacement", 0.6),
        std::make_pair<std::string, scal>("Ground029_4K_Displacement", 0.7),
        std::make_pair<std::string, scal>("Bark010_4K_Displacement", 0.7),
        std::make_pair<std::string, scal>("Concrete018_4K_Displacement", 0.8),
        std::make_pair<std::string, scal>("PaintedPlaster006_4K_Displacement", 0.8),
        std::make_pair<std::string, scal>("PaintedPlaster018_4K_Displacement", 0.8),
        std::make_pair<std::string, scal>("PavingStones025_4K_Displacement", 0.8),
        std::make_pair<std::string, scal>("PavingStones040_4K_Displacement", 0.8),
        std::make_pair<std::string, scal>("Rocks022_4K_Displacement", 0.8),
        std::make_pair<std::string, scal>("Concrete018_4K_Displacement", 0.8),
        std::make_pair<std::string, scal>("Bricks084_4K_Displacement", 0.8),
        std::make_pair<std::string, scal>("PavingStones121_4K_Displacement", 0.8),
        std::make_pair<std::string, scal>("Bricks028_4K_Displacement", 0.9),
        std::make_pair<std::string, scal>("Gravel031_4K_Displacement", 0.9),
        std::make_pair<std::string, scal>("Gravel033_4K_Displacement", 0.9),
        std::make_pair<std::string, scal>("Gravel034_4K_Displacement", 0.9),
        std::make_pair<std::string, scal>("Ground019_4K_Displacement", 0.9),
        std::make_pair<std::string, scal>("Ground020_4K_Displacement", 0.9),
        std::make_pair<std::string, scal>("TactilePaving006_4K_Displacement", 0.9),
        std::make_pair<std::string, scal>("PavingStones070_4K_Displacement", 0.9)
    };
    auto scale_iterator = std::find_if(scales.begin(), scales.end(), [mesh](std::pair<std::string, scal> scale) { return scale.first == mesh->name; });
    if (scale_iterator != scales.end()) {
        scal rdm = ((float)rand()) / (float)RAND_MAX; // between 0 and 1
        rdm = rdm * 0.5 + 0.5; // between 0.5 and 1
        scal factor = scale_iterator->second * rdm;
        Console::print(OutLevel::TRACE, Console::timeStamp.str() + "Flattening by " + std::to_string(factor) + "...");
        mesh->scale({ 1., 1., factor });
    }
    */

    Console::print(OutLevel::TRACE, Console::timeStamp.str() + "Rotation...");
    mesh->rotate(normalize(vec3sc(0, 0, 1)));

    Console::print(OutLevel::TRACE, Console::timeStamp.str() + "Scale to 64x64...");
    scal size = std::max({ mesh->bounds.upper.x - mesh->bounds.lower.x, mesh->bounds.upper.y - mesh->bounds.lower.y });
    scal halfSize = size / 2.;
    mesh->scale(64. / halfSize);

    Console::print(OutLevel::TRACE, Console::timeStamp.str() + "Shift to center and min at 0...");
    mesh->translate({ -(mesh->bounds.lower.x + mesh->bounds.upper.x) / (scal)2.,
                      -(mesh->bounds.lower.y + mesh->bounds.upper.y) / (scal)2.,
                       -mesh->bounds.lower.z });

    std::stringstream ss; ss << *mesh;
    Console::print(OutLevel::TRACE, *mesh);

    return mesh;
}


TriangleMesh* createFlatGrid()
{
    TriangleMesh* mesh = new TriangleMesh();
    gdt::vec3sc center = { 0, 0, 0 };
    scal size = 20;
    gdt::vec3sc start = { center.x - size / (scal)2., center.y - size / (scal)2., center.z };
    int n = 30;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            scal x = start.x + i * size / (scal)(n - 1);
            scal y = start.y + j * size / (scal)(n - 1);
            scal z = static_cast <scal> (rand()) / static_cast <scal> (RAND_MAX);
            z = (float)((int)(z * 100 + .5))/100;
            mesh->vertex.push_back({ x, y, z });
            mesh->bounds.extend(mesh->vertex[mesh->vertex.size() - 1]);
            if (i < n - 1 && j < n - 1) {
                mesh->index.push_back({ i * n + j, i * n + j + 1, (i + 1) * n + j });
                mesh->index.push_back({ i * n + j + 1, (i + 1) * n + j, (i + 1) * n + j + 1 });
            }
        }
    }
    mesh->meso_normal = mesh->computeMesoNormal(0);
    mesh->macroArea = mesh->bounds.span().x * mesh->bounds.span().y;
    mesh->surfaceArea = 0;
    for (int faceID = 0; faceID < mesh->index.size(); ++faceID) {
        mesh->triangle_normal.push_back(triangleNormal(mesh, mesh->index[faceID]));
        mesh->area.push_back(computeArea(mesh, mesh->index[faceID]));
        mesh->surfaceArea += mesh->area[mesh->area.size() - 1];
    }

    scal halfSize = size / 2.;
    mesh->scale(64. / halfSize);

    return mesh;
}


std::ostream& operator<<(std::ostream& output, const TriangleMesh& mesh) {
    output << Console::timePad << "Triangle mesh: " << mesh.name << std::endl;
    output << Console::timePad << Console::indent << mesh.vertex.size() << " vertices / " << mesh.index.size() << " faces" << std::endl;
    output << Console::timePad << Console::indent << "Bounds: " << mesh.bounds << std::endl;
    output << Console::timePad << Console::indent << "Macrourface area: " << mesh.macroArea << " / Microsurface area: " << mesh.surfaceArea << std::endl;
    output << Console::timePad << Console::indent << "Meso normal: " << mesh.meso_normal << std::endl;
    return output;
}