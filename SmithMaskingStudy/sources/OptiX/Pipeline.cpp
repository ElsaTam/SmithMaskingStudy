#include "OptiX/Pipeline.h"

#include <fstream>
#include <stdlib.h>
#include "OptiX/Context.h"
#include "utils/console.h"
#include "utils/params.h"


extern "C" char embedded_ptx_code[];

static bool readFile(std::string& str, const char* path)
{
    std::ifstream fp(path);
    if (fp.good())
    {
        std::stringstream content;
        content << fp.rdbuf();
        str = content.str();
        return true;
    }
    return false;
}


Pipeline::Pipeline(const char* ptx_path_)
    :
    maxTraceDepth(2),
    numPayloadValues(2),
    numAttributeValues(2)
{
    createOptions();
    createModule(ptx_path_);
    createPG();
    linkPipeline();
}

void Pipeline::createOptions()
{
    createModuleCompileOptions();
    createProgramGroupOptions();
    createPipelineOptions();
    createPipelineLinkOptions();
}

void Pipeline::createModuleCompileOptions()
{
    moduleCompileOptions = {};
    moduleCompileOptions.maxRegisterCount = 50;
    moduleCompileOptions.optLevel = OPTIX_COMPILE_OPTIMIZATION_DEFAULT;
    moduleCompileOptions.debugLevel = OPTIX_COMPILE_DEBUG_LEVEL_NONE;
    moduleCompileOptions.boundValues = NULL;
    moduleCompileOptions.numBoundValues = 0;
    moduleCompileOptions.numPayloadTypes = 0;
}

void Pipeline::createProgramGroupOptions()
{
    programGroupOptions = {}; // Initialize to zeros
}

void Pipeline::createPipelineOptions()
{
    pipelineCompileOptions = {};
    pipelineCompileOptions.traversableGraphFlags = OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_GAS;
    pipelineCompileOptions.usesMotionBlur = false;
    pipelineCompileOptions.numPayloadValues = numPayloadValues;
    pipelineCompileOptions.numAttributeValues = numAttributeValues;
    pipelineCompileOptions.exceptionFlags = OPTIX_EXCEPTION_FLAG_STACK_OVERFLOW;
    pipelineCompileOptions.pipelineLaunchParamsVariableName = "optixLaunchParams";
}

void Pipeline::createPipelineLinkOptions()
{
    pipelineLinkOptions.maxTraceDepth = maxTraceDepth;
}


/**
Pipeline::createModule
-------------------

PTX from file is read and compiled into the module

**/

void Pipeline::createModule(const char* ptx_path) 
{
    std::string ptxCode;
    if (!readFile(ptxCode, ptx_path)) {
        Console::err << "Error reading the file " << ptx_path << std::endl;
        exit(2);
    }

    if (Parameters::get()->currentParams()->outLevel >= OutLevel::INFO) {
        Console::out
            << Console::timePad << Console::indent << "ptx_path " << ptx_path << std::endl
            << Console::timePad << Console::indent << "ptx size " << ptxCode.size() << std::endl;
    }

    char log[2048]; // For error reporting from OptiX creation functions
    size_t sizeof_log = sizeof(log);

    OPTIX_CHECK_LOG(optixModuleCreateFromPTX(
        Context::optixContext,
        &moduleCompileOptions,
        &pipelineCompileOptions,
        ptxCode.c_str(),
        ptxCode.size(),
        log, &sizeof_log,
        &module
    ));
}



void Pipeline::createPG()
{
    createRaygenPG("globalVisibility");

    CudaMissPrograms cudaMissPrograms; // Miss (MS)
    cudaMissPrograms[RayType::RADIANCE_RAY_TYPE] = "radiance";
    cudaMissPrograms[RayType::AMBIENT_OCCLUSION_RAY_TYPE] = "visibility";
    cudaMissPrograms[RayType::VISIBILITY_RAY_TYPE] = "visibility";
    createMissPG(cudaMissPrograms);

    CudaHitPrograms cudaHitPrograms; // Intersection (IS), Closest Hit (CH), Any Hit (AH)
    cudaHitPrograms[RayType::RADIANCE_RAY_TYPE] = { nullptr, "radiance", nullptr };
    cudaHitPrograms[RayType::AMBIENT_OCCLUSION_RAY_TYPE] = { nullptr, "ambient_occlusion", nullptr };
    cudaHitPrograms[RayType::VISIBILITY_RAY_TYPE] = { nullptr, nullptr, "visibility" };
    createHitgroupPG(cudaHitPrograms);
}

/**
Pipeline::createRaygenPG
---------------------

Creates member raygen_pg

**/

void Pipeline::createRaygenPG(const char* rg)
{
    raygenPGs.resize(1);

    std::string rg_ = "__raygen__";
    rg_ += rg;

    OptixProgramGroupDesc pgDesc = {};
    pgDesc.kind = OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
    pgDesc.raygen.module = module;
    pgDesc.raygen.entryFunctionName = rg_.c_str();

    size_t sizeof_log = 0;
    char log[2048];

    OPTIX_CHECK_LOG(optixProgramGroupCreate(
        Context::optixContext,
        &pgDesc,
        1,
        &programGroupOptions,
        log, &sizeof_log,
        &raygenPGs[0]
    ));
}

/**
Pipeline::createMissPG
---------------------

Creates member miss_pg

**/

void Pipeline::createMissPG(CudaMissPrograms cudaPrograms)
{
    missPGs.resize(cudaPrograms.size());

    OptixProgramGroupDesc pgDesc = {};
    pgDesc.kind = OPTIX_PROGRAM_GROUP_KIND_MISS;
    pgDesc.miss.module = module;

    for (auto& ms : cudaPrograms)
    {
        RayType rayType = ms.first;

        std::string ms_ = "__miss__";
        ms_ += ms.second;

        pgDesc.miss.entryFunctionName = ms_.c_str();

        size_t sizeof_log = 0;
        char log[2048];

        OPTIX_CHECK_LOG(optixProgramGroupCreate(
            Context::optixContext,
            &pgDesc,
            1,
            &programGroupOptions,
            log, &sizeof_log,
            &missPGs[rayType]
        ));
    }
}

/**
Pipeline::createHitgroupPG
---------------------

Creates member hitgroup_pg

**/

void Pipeline::createHitgroupPG(CudaHitPrograms cudaPrograms)
{
    hitgroupPGs.resize(cudaPrograms.size());
    for (auto& hg : cudaPrograms)
    {
        OptixProgramGroupDesc pgDesc = {};
        pgDesc.kind = OPTIX_PROGRAM_GROUP_KIND_HITGROUP;

        RayType rayType = hg.first;
        const char* is = std::get<0>(hg.second);
        const char* ch = std::get<1>(hg.second);
        const char* ah = std::get<2>(hg.second);

        if (is)
        {
            char is_[100]; strcpy_s(is_, _countof(is_), "__intersection__"); strcat_s(is_, _countof(is_), is);
            pgDesc.hitgroup.moduleIS = module;
            pgDesc.hitgroup.entryFunctionNameIS = is_;
        }
        if (ch)
        {
            char ch_[100]; strcpy_s(ch_, _countof(ch_), "__closesthit__"); strcat_s(ch_, _countof(ch_), ch);
            pgDesc.hitgroup.moduleCH = module;
            pgDesc.hitgroup.entryFunctionNameCH = ch_;
        }
        if (ah)
        {
            char ah_[100]; strcpy_s(ah_, _countof(ah_), "__anyhit__"); strcat_s(ah_, _countof(ah_), ah);
            pgDesc.hitgroup.moduleAH = module;
            pgDesc.hitgroup.entryFunctionNameAH = ah_;
        }

        char log[2048];
        size_t sizeof_log = sizeof(log);

        OPTIX_CHECK_LOG(optixProgramGroupCreate(
            Context::optixContext,
            &pgDesc,
            1,
            &programGroupOptions,
            log, &sizeof_log,
            &hitgroupPGs[rayType]
        ));
    }
}

/**
Pipeline::linkPipeline
-------------------

Create pipeline from the program_groups

**/

void Pipeline::linkPipeline()
{
    std::vector<OptixProgramGroup> programGroups;
    for (auto pg : raygenPGs)
        programGroups.push_back(pg);
    for (auto pg : missPGs)
        programGroups.push_back(pg);
    for (auto pg : hitgroupPGs)
        programGroups.push_back(pg);

    size_t sizeof_log = 0;
    char log[2048];

    OPTIX_CHECK_LOG(optixPipelineCreate(
        Context::optixContext,
        &pipelineCompileOptions,
        &pipelineLinkOptions,
        programGroups.data(),
        (int)programGroups.size(),
        log, &sizeof_log,
        &pipeline
    ));
}