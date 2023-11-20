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

#pragma once

#include <queue>
#include <set>
#include <vector>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/string.hpp>

#include "gdt/math/AffineSpace.h"
  
struct TriangleMesh {
    ~TriangleMesh(){ }
    
    void scale(gdt::vec3sc factor);
    void translate(gdt::vec3sc vector);
    void rotate(gdt::vec3sc newNormal);
    int subdivisions() const;

    void createPatches();

    /**
    * @param faceID id of the triangle in the middle of the patch
    * @param d      distance to the border of the patch (a percentage of the surface size)
    */
    gdt::vec3sc neighborhoodMesoNormal(int faceID, float d) const;
    gdt::vec3sc neighborhoodMesoNormal(int faceID) const;

    gdt::vec3sc computeMesoNormal(int maxSamples) const;


    void flatSetContainingFace(int faceID, gdt::vec3sc N, std::queue<int>& faceIDQueue) const;
    std::vector<std::vector<int>> flatSets(gdt::vec3sc N) const;

    std::vector<std::set<int>> heightSeparation(int K, std::vector<scal>& centroids) const;

    friend std::ostream& operator<<(std::ostream& output, const TriangleMesh& mesh);


    std::vector<gdt::vec3sc> vertex;
    std::vector<gdt::vec3sc> vertex_normal; // one normal per vertex
    std::vector<gdt::vec3sc> triangle_normal; // one normal per triangle
    std::vector<gdt::vec3i> index;
    std::vector<scal> area; // one area per triangle
    scal surfaceArea = 0;
    scal macroArea = 0;
    scal maxArea = 0; // area of the largest triangle
    gdt::vec3sc meso_normal;

    //! bounding box of all vertices in the mesh
    gdt::box3sc bounds;
    std::vector< std::pair<gdt::box3sc, gdt::vec3sc> > patches; // one bounding box associated to one normal
    std::string name;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & vertex;
        ar & vertex_normal;
        ar & triangle_normal;
        ar & index;
        ar & area;
        ar & surfaceArea;
        ar & maxArea;
        ar & meso_normal;
        ar & bounds;
        ar & patches;
        ar & name;
    }
};

TriangleMesh* createMesh(const std::string& objPath);
TriangleMesh* loadOBJ(const std::string& objPath);
TriangleMesh* createSubMesh(const TriangleMesh* _mesh, std::set<int> faceIDs);
gdt::vec3sc fittingPlaneNormal(const std::vector<gdt::vec3sc>& vertices);
TriangleMesh* createFlatGrid();