#ifndef BLADE_BASE_HH
#define BLADE_BASE_HH

#include "blade/types.hh"
#include "blade/macros.hh"
#include "blade/logger.hh"
#include "blade/pipeline.hh"

namespace Blade {

inline const Result SetCudaDevice(int device_id) {
    BL_CUDA_CHECK(cudaSetDevice(device_id), [&]{
       BL_FATAL("Failed to set device: {}", err);
    });
    return Result::SUCCESS;
}

}  // namespace Blade

#endif
