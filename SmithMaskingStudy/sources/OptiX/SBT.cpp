#include "OptiX/SBT.h"
#include "OptiX/Binding.h"
#include "utils/console.h"


/**
SBT
====

SBT needs PIP as the packing of SBT record headers requires
access to their corresponding program groups (PGs).
This is one aspect of establishing the connection between the
PGs and their data.

**/

SBT::SBT(const Pipeline* pip_)
    :
    pipeline(pip_)
{
    buildRaygenRecords();
    buildMissRecords();
}



void SBT::buildRaygenRecords()
{
    std::vector<RaygenRecord> raygenRecords;

    RaygenRecord rec;
    OPTIX_CHECK(optixSbtRecordPackHeader(pipeline->getRaygenPGs()[0], &rec));
    rec.data.data = nullptr; /* for now ... */
    raygenRecords.push_back(rec);

    if (raygenRecordsBuffer.d_ptr != nullptr)
        raygenRecordsBuffer.freeBuffer();
    raygenRecordsBuffer.alloc_and_upload(raygenRecords);

    sbt.raygenRecord = raygenRecordsBuffer.d_pointer();
}

void SBT::buildMissRecords()
{
    std::vector<MissRecord> missRecords;

    for (int i = 0; i < pipeline->getMissPGs().size(); i++) {
        MissRecord rec;
        OPTIX_CHECK(optixSbtRecordPackHeader(pipeline->getMissPGs()[i], &rec));
        rec.data.data = nullptr; /* for now ... */
        missRecords.push_back(rec);
    }

    if (missRecordsBuffer.d_ptr != nullptr)
        missRecordsBuffer.freeBuffer();
    missRecordsBuffer.alloc_and_upload(missRecords);

    sbt.missRecordBase = missRecordsBuffer.d_pointer();
    sbt.missRecordStrideInBytes = sizeof(MissRecord);
    sbt.missRecordCount = (int)missRecords.size();
}



void SBT::updateAS(const AccelerationStructure& AS)
{
    buildHitgroupRecords(AS);
}

void SBT::buildHitgroupRecords(const AccelerationStructure& AS)
{
    std::vector<HitgroupRecord> hitgroupRecords;

    for (int rayID = 0; rayID < RAY_TYPE_COUNT; rayID++) {
        HitgroupRecord rec;
        OPTIX_CHECK(optixSbtRecordPackHeader(pipeline->getHitgroupPGs()[rayID], &rec));
        rec.data.triangleMesh.vertex = (vec3sc*)AS.getVertexBuffer().d_pointer();
        rec.data.triangleMesh.normal = (vec3sc*)AS.getNormalBuffer().d_pointer();
        rec.data.triangleMesh.index = (vec3i*)AS.getIndexBuffer().d_pointer();
        rec.data.triangleMesh.bounds = AS.getBounds();
        Console::info << rec.data.triangleMesh.bounds.span() << std::endl;
        hitgroupRecords.push_back(rec);
    }

    if (hitgroupRecordsBuffer.d_ptr != nullptr)
        hitgroupRecordsBuffer.freeBuffer();
    hitgroupRecordsBuffer.alloc_and_upload(hitgroupRecords);

    sbt.hitgroupRecordBase = hitgroupRecordsBuffer.d_pointer();
    sbt.hitgroupRecordStrideInBytes = sizeof(HitgroupRecord);
    sbt.hitgroupRecordCount = (int)hitgroupRecords.size();
}