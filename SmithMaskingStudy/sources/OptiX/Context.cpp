#include "OptiX/Context.h"
#include "utils/console.h"
#include "utils/params.h"

OptixDeviceContext Context::optixContext = nullptr;

void Context::context_log_cb(unsigned int level, const char* tag, const char* message, void* /*cbdata */)  // static 
{
    // callback levels: OutLevel
    // 0: NO_OUTPUT
    // 1: ERR
    // 2: WARNING
    // 3: INFO
    // 4: SUCCESS
    // 5: NORMAL
    // 6: TRACE

    std::stringstream ss;
    ss << Console::timePad;
    ss << "[" << std::setw(2) << level << "][" << std::setw(12) << tag << "]: " << message;
    Console::printOptix(OutLevel(level), ss.str());
}



Context::Context()
{
    if (!optixContext)
    {

        Console::print(OutLevel::TRACE, Console::timeStamp.str() + "Creating OptiX context...");

        // -------------------------------------------------------
        // init optix
        // -------------------------------------------------------

        cudaFree(0);
        int numDevices;
        cudaGetDeviceCount(&numDevices);
        if (numDevices == 0)
            throw std::runtime_error("Optix: no CUDA capable devices found!");
        OPTIX_CHECK(optixInit());
        Console::print(OutLevel::SUCCESS, Console::timePad + "Optix successgully initialized.");


        // -------------------------------------------------------
        // create context
        // -------------------------------------------------------

        // do everything on one device
        const int deviceID = 0;
        CUDA_CHECK(SetDevice(deviceID));

        // check which device we are working on
        cudaDeviceProp deviceProps;
        cudaGetDeviceProperties(&deviceProps, deviceID);
        Console::print(OutLevel::NORMAL, Console::timePad + "Optix running on device: " + deviceProps.name);

        // Get the current cuda context
        CUresult cuRes = cuCtxGetCurrent(&cudaContext);
        if (cuRes != CUDA_SUCCESS) {
            std::stringstream ss_cuRes;
            ss_cuRes << cuRes;
            Console::print(OutLevel::ERR, "Error querying current context: error code " + ss_cuRes.str());
        }

        // Set the callback
        OptixDeviceContextOptions options = {};
        options.logCallbackFunction = &Context::context_log_cb;
        options.logCallbackLevel = (int) OutLevel::NORMAL;

        // Create the context
        OPTIX_CHECK(optixDeviceContextCreate(cudaContext, &options, &optixContext));
    }
}