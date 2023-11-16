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

#include <optix_device.h>

#include <OptiX/LaunchParams.h>
#include <gdt/random/random.h>
#include <gdt/random/sampling.h>

/*! launch parameters in constant memory, filled in by optix upon
      optixLaunch (this gets filled in from the buffer we pass to
      optixLaunch) */
extern "C" __constant__ LaunchParams optixLaunchParams;

/*! per-ray data now captures random number generator, so programs
      can access RNG state */
struct PRD {
    Random<16> random;
    unsigned int seed;
    vec3sc color;
    bool backface;
    bool visible;
    bool onSurface;
};

static __forceinline__ __device__
void* unpackPointer(uint32_t i0, uint32_t i1)
{
    const uint64_t uptr = static_cast<uint64_t>(i0) << 32 | i1;
    void* ptr = reinterpret_cast<void*>(uptr);
    return ptr;
}

static __forceinline__ __device__
void  packPointer(void* ptr, uint32_t& i0, uint32_t& i1)
{
    const uint64_t uptr = reinterpret_cast<uint64_t>(ptr);
    i0 = uptr >> 32;
    i1 = uptr & 0x00000000ffffffff;
}

template<typename T>
static __forceinline__ __device__ T* getPRD()
{
    const uint32_t u0 = optixGetPayload_0();
    const uint32_t u1 = optixGetPayload_1();
    return reinterpret_cast<T*>(unpackPointer(u0, u1));
}

__forceinline__ __device__ vec3sc AABBIntersection(const box3sc& aabb, const vec3sc& O, const vec3sc& D)
{
    vec3sc v_tmin = vec3sc(
        (aabb.lower.x - O.x) / D.x,
        (aabb.lower.y - O.y) / D.y,
        (aabb.lower.z - O.z) / D.z
    );
    vec3sc v_tmax = vec3sc(
        (aabb.upper.x - O.x) / D.x,
        (aabb.upper.y - O.y) / D.y,
        (aabb.upper.z - O.z) / D.z
    );
    for (int i = 0; i < 3; ++i) {
        if (v_tmin[i] > v_tmax[i]) {
            scal tmp = v_tmin[i];
            v_tmin[i] = v_tmax[i];
            v_tmax[i] = tmp;
        }
    }
    scal t = min(v_tmax[0], min(v_tmax[1], v_tmax[2]));
    return O + t * D;
}

//------------------------------------------------------------------------------
// closest hit and anyhit programs for radiance-type rays.
//
// Note eventually we will have to create one pair of those for each
// ray type and each geometry type we want to render; but this
// simple example doesn't use any actual geometries yet, so we only
// create a single, dummy, set of them (we do have to have at least
// one group of them to set up the SBT)
//------------------------------------------------------------------------------

// closest hit program for a vertical ray launched somewhere above the surface.
// we want to find the point on the surface juste under the ray origin.
// we need to take care of the side effect and discard unreliable points.
extern "C" __global__ void __closesthit__radiance()
{
    // get the per-ray data structure reference :
    PRD& prd = *getPRD<PRD>();

    // get the mesh data :
    const TriangleMeshSBTData& sbtData
        = *(const TriangleMeshSBTData*)optixGetSbtDataPointer();

    // get the triangle data :
    const int   primID = optixGetPrimitiveIndex(); // faceID
    const vec3i index = sbtData.index[primID];
    const vec3sc& A = sbtData.vertex[index.x];
    const vec3sc& B = sbtData.vertex[index.y];
    const vec3sc& C = sbtData.vertex[index.z];
    const scal u = optixGetTriangleBarycentrics().x;
    const scal v = optixGetTriangleBarycentrics().y;

    // check if face belongs to a set
    if (optixLaunchParams.sets.numberSets > 0) {
        prd.color = optixLaunchParams.sets.colorsBySet[optixLaunchParams.sets.colorsId[primID]];
    }

    // compute the geometric normal :
    vec3sc N = (optixLaunchParams.visibility.useSmooth && sbtData.normal) ?
        ((1.f - u - v) * sbtData.normal[index.x]
            + u * sbtData.normal[index.y]
            + v * sbtData.normal[index.z])
        : normalize(cross(B - A, C - A));
    if (N.z < 0.f) N = -N;

    vec3sc rayDir;
    unsigned int flags;
    unsigned int rayType;
    if (   optixLaunchParams.camera.programType == ProgramType::G1
        || optixLaunchParams.camera.programType == ProgramType::GAF)
    {
        rayDir = optixLaunchParams.visibility.directionOut;
        flags = OPTIX_RAY_FLAG_TERMINATE_ON_FIRST_HIT | OPTIX_RAY_FLAG_DISABLE_CLOSESTHIT;
        rayType = VISIBILITY_RAY_TYPE;
    }
    else if (optixLaunchParams.camera.programType == ProgramType::AMBIENT_OCCLUSION)
    {
        scal r1 = prd.random.rnd(prd.seed); // [0, 1]
        scal r2 = prd.random.rnd(prd.seed); // [0, 1]
        vec3sc verticalDir = uniformHemisphereSampling(r1, r2);
        rayDir = rotateAlongNormal(verticalDir, N);

        flags = OPTIX_RAY_FLAG_NONE;
        rayType = AMBIENT_OCCLUSION_RAY_TYPE;
    }
    else
    {
        return;
    }

    // find the 3D point on the surface.
    const vec3sc surfPos = (1.f - u - v) * A + u * B + v * C;

    // If we are too close to the edge of the surface, discard the point for the minimal visibility
    // (if optixLaunchParams.sideEffect.borderPercentage > 0)
    if (sbtData.bounds.closest_distance(surfPos).x < optixLaunchParams.sideEffect.borderPercentage * sbtData.bounds.span().x / 2.f
        || sbtData.bounds.closest_distance(surfPos).y < optixLaunchParams.sideEffect.borderPercentage * sbtData.bounds.span().y / 2.f)
    {
        prd.onSurface = false;
        return;
    }
    prd.onSurface = true;

    if (dot(N, rayDir) <= 0) {
        prd.backface = true;
        return;
    }
    prd.backface = false;

    // pack the visibility boolean :
    uint32_t u0, u1;
    packPointer(&(prd.visible), u0, u1);

    // launch a ray in the visibility direction (OUTgoing) :
    optixTrace(optixLaunchParams.traversable,
        surfPos + (scal)1e-3 * N,
        rayDir,
        (scal)0.0, // tmin
        optixLaunchParams.visibility.tMax,  // tmax
        (scal)0.0,  // rayTime
        OptixVisibilityMask(255),
        flags,
        rayType,            // SBT offset
        RAY_TYPE_COUNT,     // SBT stride
        rayType,            // missSBTIndex 
        u0, u1);

    // If the point is not visible, then there is no possible mistake: the surface is hiding it.
    if (optixLaunchParams.sideEffect.directional && prd.visible)
    {
        vec3sc hit = AABBIntersection(sbtData.bounds, surfPos, -optixLaunchParams.visibility.directionOut);
        if (hit.z < sbtData.bounds.upper.z && hit.z > sbtData.bounds.lower.z) {
            prd.onSurface = false;
            return;
        }
    }

    if (prd.visible && optixLaunchParams.camera.programType == ProgramType::GAF) {
        rayDir = optixLaunchParams.visibility.directionIn;
        if (dot(N, rayDir) <= 0) {
            prd.visible = false;
            return;
        }
        optixTrace(optixLaunchParams.traversable,
            surfPos + (scal)1e-3 * N,
            rayDir,
            (scal)0.0, // tmin
            optixLaunchParams.visibility.tMax,  // tmax
            (scal)0.0,  // rayTime
            OptixVisibilityMask(255),
            flags,
            VISIBILITY_RAY_TYPE, // SBT offset
            RAY_TYPE_COUNT,      // SBT stride
            VISIBILITY_RAY_TYPE, // missSBTIndex 
            u0, u1);

        if (optixLaunchParams.sideEffect.directional && prd.visible)
        {
            vec3sc hit = AABBIntersection(sbtData.bounds, surfPos, -optixLaunchParams.visibility.directionIn);
            if (hit.z < sbtData.bounds.upper.z && hit.z > sbtData.bounds.lower.z) {
                prd.onSurface = false;
                return;
            }
        }
    }
}

// closest hit program vor a directional ray launched from the surface,
// with a specified tMax.
extern "C" __global__ void __closesthit__ambient_occlusion()
{
    // we hitted something, so the point is not visible
    bool& visibility = *getPRD<bool>();
    visibility = false;
}


//------------------------------------------------------------------------------
// any hit programs.
//------------------------------------------------------------------------------

extern "C" __global__ void __anyhit__visibility() {
    const TriangleMeshSBTData& sbtData = *(const TriangleMeshSBTData*)optixGetSbtDataPointer();
    const vec3i index = sbtData.index[optixGetPrimitiveIndex()];
    const scal u = optixGetTriangleBarycentrics().x;
    const scal v = optixGetTriangleBarycentrics().y;
    const vec3sc hit = (1.f - u - v) * sbtData.vertex[index.x] + u * sbtData.vertex[index.y] + v * sbtData.vertex[index.z];

    bool& visibility = *getPRD<bool>();

    // If we are too close to the edge of the surface, consider no hit (visible)
    const scal epsilonX = 0.001 * sbtData.bounds.span().x;
    const scal epsilonY = 0.001 * sbtData.bounds.span().y;
    if (hit.x <= sbtData.bounds.lower.x + epsilonX || hit.x >= sbtData.bounds.upper.x - epsilonX
        || hit.y <= sbtData.bounds.lower.y + epsilonY || hit.y >= sbtData.bounds.upper.y - epsilonY)
    {
        visibility = true;
        return;
    }

    // we hitted something, so the point is not visible
    visibility = false;
}


//------------------------------------------------------------------------------
// miss programs.
// ------------------------------------------------------------------------------

extern "C" __global__ void __miss__radiance()
{
    PRD& prd = *getPRD<PRD>();
    prd.onSurface = false;
}

extern "C" __global__ void __miss__visibility()
{
    // we didn't hit anything, so the point is visible
    bool& visibility = *getPRD<bool>();
    visibility = true;
}



//------------------------------------------------------------------------------
// ray generation programs.
//------------------------------------------------------------------------------

extern "C" __global__ void __raygen__globalVisibility()
{
    // compute a test pattern based on pixel ID
    const scal ix = optixGetLaunchIndex().x;
    const scal iy = optixGetLaunchIndex().y;
    const int accumID = optixLaunchParams.frame.accumID;
    const auto& camera = optixLaunchParams.camera;

    // our per-ray data for this example. what we initialize it to
    // won't matter, since this value will be overwritten by either
    // the miss or hit program, anyway
    PRD prd;
    prd.seed = prd.random.init(ix + accumID * optixLaunchParams.frame.size.x,
        iy + accumID * optixLaunchParams.frame.size.y);
    prd.color = { 0.f, 0.f, 0.f };
    prd.backface = false;
    prd.visible = false;
    prd.onSurface = false;
    vec3sc pixelColorPRD{ 0.f, 0.f, 0.f };
    uint32_t validRays = 0;
    uint32_t visibleRays = 0;

    // the values we store the PRD pointer in:
    uint32_t u0, u1;
    packPointer(&prd, u0, u1);

    int numPixelSamples = optixLaunchParams.camera.nPixelSamples;

    for (int sampleID = 0; sampleID < numPixelSamples; sampleID++) {
        // normalized screen plane position, in [0,1]^2
        const float screen_u1 = prd.random.rnd(prd.seed);
        const float screen_u2 = prd.random.rnd(prd.seed);
        const vec2sc screen = vec2sc(ix + screen_u1, iy + screen_u2) / vec2sc(optixLaunchParams.frame.size);

        // generate ray origin (orthographic camera)
        vec3sc rayOrigin = camera.position
            + (screen.x - (scal)0.5) * camera.horizontal
            + (screen.y - (scal)0.5) * camera.vertical;

        optixTrace(optixLaunchParams.traversable,
            rayOrigin,
            optixLaunchParams.camera.direction, // ray direction
            (scal)0,    // tmin
            (scal)1e20,  // tmax
            (scal)0.0,   // rayTime
            OptixVisibilityMask(255),
            OPTIX_RAY_FLAG_DISABLE_ANYHIT,  // OPTIX_RAY_FLAG_NONE,
            RADIANCE_RAY_TYPE,            // SBT offset
            RAY_TYPE_COUNT,               // SBT stride
            RADIANCE_RAY_TYPE,            // missSBTIndex 
            u0, u1);

        if (optixLaunchParams.sets.numberSets > 0) {
            pixelColorPRD += prd.color;
        }
        else if (prd.onSurface && !prd.backface) {
            validRays++;
            if (prd.visible) {
                pixelColorPRD += vec3sc(0.f, 1.f, 0.f);
                visibleRays++;
            }
            else {
                pixelColorPRD += vec3sc(1.f, 0.f, 0.f);
            }
        }
        else {
            pixelColorPRD += vec3sc(0.f, 0.f, 1.f);
        }
    }

    pixelColorPRD /= numPixelSamples;
    int r, g, b;
    if (optixLaunchParams.camera.programType == ProgramType::AMBIENT_OCCLUSION) {
        r = g = b = int(255.99f * (scal)visibleRays / (scal)validRays);
    }
    else {
        r = int(255.99f * pixelColorPRD.x);
        g = int(255.99f * pixelColorPRD.y);
        b = int(255.99f * pixelColorPRD.z);
    }
    const uint32_t rgba = 0xff000000
        | (r << 0) | (g << 8) | (b << 16);


    // and write to frame buffer ...
    const uint32_t fbIndex = ix + iy * optixLaunchParams.frame.size.x;

    optixLaunchParams.frame.colorBuffer[fbIndex] = rgba;
    if (validRays == 0) {
        optixLaunchParams.frame.visibilityBuffer[fbIndex] = -1.f;
    }
    else {
        optixLaunchParams.frame.visibilityBuffer[fbIndex] = (scal)visibleRays / (scal)validRays;
    }
}