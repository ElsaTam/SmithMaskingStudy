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

#include <vector>
#include "gdt/math/box.h"
#include "optix7.h" // for OptixTraversableHandle

using namespace gdt;

enum RayType { RADIANCE_RAY_TYPE = 0, AMBIENT_OCCLUSION_RAY_TYPE, VISIBILITY_RAY_TYPE, RAY_TYPE_COUNT };
enum class ProgramType { G1 = 0, GAF, AMBIENT_OCCLUSION};

struct TriangleMeshSBTData {
    vec3sc* vertex;
    vec3sc* normal;
    vec3i* index;
    box3sc  bounds;
};

struct LaunchParams
{
    struct {
        uint32_t* colorBuffer;
        scal*    visibilityBuffer;
        vec2i     size;
        int       accumID{ 0 };
    } frame;

    struct {
        float cosFovy;
        vec3sc position;
        vec3sc direction;
        vec3sc horizontal;
        vec3sc vertical;
        int nPixelSamples;
        ProgramType programType;
    } camera;

    struct {
        vec3sc directionIn;
        vec3sc directionOut;
        size_t* targetedTriangles;
        size_t numberOfTargetedTriangles{ 0 };
        scal tMax;
        bool useSmooth;
    } visibility;

    struct {
        vec3sc* colorsBySet;
        int* colorsId;
        size_t numberSets{ 0 };
    } sets;

    struct {
        float borderPercentage;
        bool directional;
    } sideEffect;

    OptixTraversableHandle traversable;
};
