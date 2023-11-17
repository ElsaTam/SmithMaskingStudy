#pragma once

#include "OptiX/optix7.h"

struct Context
{
    CUcontext cudaContext;
    static OptixDeviceContext optixContext;

    static void context_log_cb(unsigned int level, const char* tag, const char* message, void* /*cbdata */);

    Context();
};

