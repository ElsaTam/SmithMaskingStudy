#pragma once

#include "OptiX/CUDABuffer.h"
#include "shapes/TriangleMesh.h"

class AccelerationStructure
{
private:
	OptixTraversableHandle asHandle;

    OptixAccelBuildOptions accelOptions;
    OptixAccelBufferSizes blasBufferSizes;

    uint32_t triangleInputFlags;

    gdt::box3sc bounds;
    CUDABuffer vertexBuffer;
    CUDABuffer normalBuffer;
    CUDABuffer indexBuffer;

    //! buffer that keeps the (final, compacted) accel structure
    CUDABuffer asBuffer;

public:
	AccelerationStructure(const TriangleMesh* mesh);

    const OptixTraversableHandle& getHandle() const { return asHandle; }
    const gdt::box3sc&     getBounds()       const { return bounds; }
    const CUDABuffer& getVertexBuffer() const { return vertexBuffer; }
    const CUDABuffer& getNormalBuffer() const { return normalBuffer; }
    const CUDABuffer& getIndexBuffer()  const { return indexBuffer; }
};