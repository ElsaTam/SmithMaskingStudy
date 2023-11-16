#include "OptiX/Context.h"
#include "utils/console.h"

OptixDeviceContext Context::optixContext = nullptr;

void Context::context_log_cb(unsigned int level, const char* tag, const char* message, void* /*cbdata */)  // static 
{
    // callback levels:
    // 0: disable
    // 1: fatal
    // 2: error
    // 3: warning
    // 4: print

    Console::OutStream out = Console::optixOut << Console::timePad;
    if (level <= 2) {
        out = Console::optixErr;
    }
    else if (level <= 3) {
        out = Console::optixWarn;
    }
    out << "[" << std::setw(2) << level << "][" << std::setw(12) << tag << "]: "
            << message << "\n";
}



Context::Context()
{
    if (!optixContext)
    {
        Console::out << Console::timeStamp << "Creating OptiX context..." << std::endl;

        // -------------------------------------------------------
        // init optix
        // -------------------------------------------------------

        cudaFree(0);
        int numDevices;
        cudaGetDeviceCount(&numDevices);
        if (numDevices == 0)
            throw std::runtime_error("Optix: no CUDA capable devices found!");
        OPTIX_CHECK(optixInit());
        Console::succ << Console::timePad << "Optix successgully initialized." << std::endl;


        // -------------------------------------------------------
        // create context
        // -------------------------------------------------------

        // do everything on one device
        const int deviceID = 0;
        CUDA_CHECK(SetDevice(deviceID));

        // check which device we are working on
        cudaDeviceProp deviceProps;
        cudaGetDeviceProperties(&deviceProps, deviceID);
        Console::light << Console::timePad << "Optix running on device: " << deviceProps.name << std::endl;

        // Get the current cuda context
        CUresult cuRes = cuCtxGetCurrent(&cudaContext);
        if (cuRes != CUDA_SUCCESS)
            Console::err << "Error querying current context: error code " << cuRes << std::endl;

        // Set the callback
        OptixDeviceContextOptions options = {};
        options.logCallbackFunction = &Context::context_log_cb;
        options.logCallbackLevel = 3;

        // Create the context
        OPTIX_CHECK(optixDeviceContextCreate(cudaContext, &options, &optixContext));
    }
}