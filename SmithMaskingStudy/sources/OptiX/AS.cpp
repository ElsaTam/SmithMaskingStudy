#include "OptiX/AS.h"
#include "OptiX/Context.h"
#include "utils/console.h"
#include "utils/params.h"

AccelerationStructure::AccelerationStructure(const TriangleMesh* mesh)
{
    // ==================================================================
    // triangle inputs
    // ==================================================================
    OptixBuildInput triangleInput;
    CUdeviceptr d_vertices;
    CUdeviceptr d_indices;
    uint32_t triangleInputFlags;

    // upload the model to the device: the builder
    const TriangleMesh& _mesh = *mesh;
    vertexBuffer.alloc_and_upload(_mesh.vertex);
    //normalBuffer.alloc_and_upload(_mesh.vertex_normal);
    indexBuffer.alloc_and_upload(_mesh.index);
    bounds = _mesh.bounds;

    Console::print(OutLevel::NORMAL, Console::timeStamp.str()
        + "-- vertexBuffer size: " + std::to_string(vertexBuffer.sizeInBytes / 1024 / 1024)
        + " MB (" + std::to_string(vertexBuffer.sizeInBytes) + " bytes)");
    Console::print(OutLevel::NORMAL, Console::timeStamp.str()
        + "-- indexBuffer size: " + std::to_string(indexBuffer.sizeInBytes / 1024 / 1024)
        + " MB (" + std::to_string(indexBuffer.sizeInBytes) + " bytes)");

    triangleInput = {};
    triangleInput.type
        = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;

    // create local variables, because we need a *pointer* to the
    // device pointers
    d_vertices = vertexBuffer.d_pointer();
    d_indices = indexBuffer.d_pointer();

    triangleInput.triangleArray.vertexFormat = OPTIX_VERTEX_FORMAT_FLOAT3;
    triangleInput.triangleArray.vertexStrideInBytes = sizeof(gdt::vec3sc);
    triangleInput.triangleArray.numVertices = (int)_mesh.vertex.size();
    triangleInput.triangleArray.vertexBuffers = &d_vertices;

    triangleInput.triangleArray.indexFormat = OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
    triangleInput.triangleArray.indexStrideInBytes = sizeof(gdt::vec3i);
    triangleInput.triangleArray.numIndexTriplets = (int)_mesh.index.size();
    triangleInput.triangleArray.indexBuffer = d_indices;

    triangleInputFlags = 0;

    // in this example we have one SBT entry, and no per-primitive
    // materials:
    triangleInput.triangleArray.flags = &triangleInputFlags;
    triangleInput.triangleArray.numSbtRecords = 1;
    triangleInput.triangleArray.sbtIndexOffsetBuffer = 0;
    triangleInput.triangleArray.sbtIndexOffsetSizeInBytes = 0;
    triangleInput.triangleArray.sbtIndexOffsetStrideInBytes = 0;

    // ==================================================================
    // BLAS setup
    // ==================================================================

    OptixAccelBuildOptions accelOptions = {};
    accelOptions.buildFlags = OPTIX_BUILD_FLAG_ALLOW_COMPACTION;
    accelOptions.motionOptions.numKeys = 1;
    accelOptions.operation = OPTIX_BUILD_OPERATION_BUILD;

    OptixAccelBufferSizes blasBufferSizes;
    OPTIX_CHECK(optixAccelComputeMemoryUsage
    (Context::optixContext,
        &accelOptions,
        &triangleInput,
        1,  // num_build_inputs
        &blasBufferSizes
    ));

    // ==================================================================
    // prepare compaction
    // ==================================================================

    CUDABuffer compactedSizeBuffer;
    compactedSizeBuffer.alloc(sizeof(uint64_t));
    Console::print(OutLevel::NORMAL, Console::timeStamp.str()
        + "-- compactedSizeBuffer size: " + std::to_string(compactedSizeBuffer.sizeInBytes / 1024 / 1024) 
        + " MB (" + std::to_string(compactedSizeBuffer.sizeInBytes) + " bytes)");

    OptixAccelEmitDesc emitDesc;
    emitDesc.type = OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
    emitDesc.result = compactedSizeBuffer.d_pointer();

    // ==================================================================
    // execute build (main stage)
    // ==================================================================

    CUDABuffer tempBuffer;
    tempBuffer.alloc(blasBufferSizes.tempSizeInBytes);
    Console::print(OutLevel::NORMAL, Console::timeStamp.str()
        + "-- tempBuffer size: " + std::to_string(tempBuffer.sizeInBytes / 1024 / 1024)
        + " MB (" + std::to_string(tempBuffer.sizeInBytes) + " bytes)");

    CUDABuffer outputBuffer;
    Console::print(OutLevel::NORMAL, Console::timeStamp.str()
        + "-- outputBuffer size: " + std::to_string(blasBufferSizes.outputSizeInBytes / 1024 / 1024)
        + " MB (" + std::to_string(blasBufferSizes.outputSizeInBytes) + " bytes)");
    outputBuffer.alloc(blasBufferSizes.outputSizeInBytes);

    OPTIX_CHECK(optixAccelBuild(Context::optixContext,
        /* stream */0,
        &accelOptions,
        &triangleInput,
        1, // one mesh
        tempBuffer.d_pointer(),
        tempBuffer.sizeInBytes,

        outputBuffer.d_pointer(),
        outputBuffer.sizeInBytes,

        &asHandle,

        &emitDesc, 1
    ));
    CUDA_SYNC_CHECK();

    // ==================================================================
    // perform compaction
    // ==================================================================
    uint64_t compactedSize;
    compactedSizeBuffer.download(&compactedSize, 1);

    asBuffer.alloc(compactedSize);
    OPTIX_CHECK(optixAccelCompact(Context::optixContext,
        /*stream:*/0,
        asHandle,
        asBuffer.d_pointer(),
        asBuffer.sizeInBytes,
        &asHandle));
    CUDA_SYNC_CHECK();
}