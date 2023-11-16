#pragma once

#include "OptiX/Pipeline.h"
#include "OptiX/AS.h"

class SBT
{
private:
    const Pipeline* pipeline;

    CUDABuffer raygenRecordsBuffer;
    CUDABuffer missRecordsBuffer;
    CUDABuffer hitgroupRecordsBuffer;

    OptixShaderBindingTable sbt = {};

public:
    SBT(const Pipeline* pip_);

    const OptixShaderBindingTable& getSBT() const { return sbt; }
    void updateAS(const AccelerationStructure& AS);

protected:
    void buildRaygenRecords();
    void buildMissRecords();
    void buildHitgroupRecords(const AccelerationStructure& AS);
};