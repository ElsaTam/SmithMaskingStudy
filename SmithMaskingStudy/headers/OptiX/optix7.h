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

#pragma once

// optix 7
#include <cuda_runtime.h>
#include <optix.h>
#include <optix_stubs.h>
#include <sstream>
#include <stdexcept>

#define CUDA_CHECK(call)							                          \
    {									                                      \
        cudaError_t rc = cuda##call;                                          \
        if (rc != cudaSuccess)                                                \
        {                                                                     \
            std::stringstream txt;                                            \
            cudaError_t err =  rc; /*cudaGetLastError();*/                    \
            txt << "CUDA Error " << cudaGetErrorName(err)                     \
                << " (" << cudaGetErrorString(err) << ")";                    \
            throw std::runtime_error(txt.str());                              \
        }                                                                     \
    }

#define CUDA_CHECK_NOEXCEPT(call)                                             \
    {									                                      \
      cuda##call;                                                             \
    }

#define CUDA_SYNC_CHECK()                                                     \
    {                                                                         \
        cudaDeviceSynchronize();                                              \
        cudaError_t error = cudaGetLastError();                               \
        if( error != cudaSuccess )                                            \
        {                                                                     \
            std::ostringstream oss;                                           \
            oss << "error (" << __FILE__ << ": line " << __LINE__ << "): "    \
                << cudaGetErrorString( error );                               \
            std::cerr << oss.str();                                           \
            throw std::runtime_error(oss.str());                              \
            exit( 2 );                                                        \
        }                                                                     \
    }

#define OPTIX_CHECK( call )                                                   \
    {                                                                         \
        OptixResult res = call;                                               \
        if( res != OPTIX_SUCCESS )                                            \
        {                                                                     \
            std::cerr << "Optix call '" << #call << "' failed with code "     \
                 << res << ": " << __FILE__ " : " << __LINE__ << ")\n";       \
            exit( 2 );                                                        \
        }                                                                     \
    }

#define OPTIX_CHECK_LOG( call )                                               \
    {                                                                         \
        OptixResult res = call;                                               \
        if( res != OPTIX_SUCCESS )                                            \
        {                                                                     \
            std::cerr << "Optix call '" << #call << "' failed with code "     \
                 << res << ": " << __FILE__ " : " << __LINE__ << ")\n"        \
                 << "Log:\n" << log                                           \
                 << ( sizeof_log > sizeof( log ) ? "<TRUNCATED>" : "" )       \
                 << "\n";                                                     \
            exit( 2 );                                                        \
        }                                                                     \
    }