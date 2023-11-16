#pragma once

#include <map>
#include "OptiX/LaunchParams.h"

class Pipeline
{
private:
    // Miss (MS)
    typedef std::map<RayType, const char*> CudaMissPrograms;
    // Intersection (IS), Closest Hit (CH), Any Hit (AH)
    typedef std::map<RayType, std::tuple<const char*, const char*, const char*>> CudaHitPrograms;

    unsigned maxTraceDepth;
    unsigned numPayloadValues;
    unsigned numAttributeValues;

    /*! @{ the pipeline we're building */
    OptixPipeline pipeline = nullptr;
    OptixPipelineCompileOptions pipelineCompileOptions = {};
    OptixPipelineLinkOptions    pipelineLinkOptions = {};

    /*! @{ the module that contains out device programs */
    OptixModule module = nullptr;
    OptixModuleCompileOptions moduleCompileOptions = {};

    OptixProgramGroupOptions programGroupOptions = {};
    std::vector<OptixProgramGroup> raygenPGs;
    std::vector<OptixProgramGroup> missPGs;
    std::vector<OptixProgramGroup> hitgroupPGs;

public:
    Pipeline(const char* ptx_path_);

    const OptixPipeline& getPipeline() const { return pipeline; }
    const std::vector<OptixProgramGroup>& getRaygenPGs() const { return raygenPGs; }
    const std::vector<OptixProgramGroup>& getMissPGs() const { return missPGs; }
    const std::vector<OptixProgramGroup>& getHitgroupPGs() const { return hitgroupPGs; }

protected:
    void createOptions();
    void createModuleCompileOptions();
    void createProgramGroupOptions();
    void createPipelineOptions();
    void createPipelineLinkOptions();

    void createModule(const char* ptx_path);

    void createPG();
    void createRaygenPG(const char* rg);
    void createMissPG(CudaMissPrograms cudaPrograms);
    void createHitgroupPG(CudaHitPrograms cudaPrograms);

    void linkPipeline();
};