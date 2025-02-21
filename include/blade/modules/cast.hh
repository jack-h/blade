#ifndef BLADE_MODULES_CAST_GENERIC_HH
#define BLADE_MODULES_CAST_GENERIC_HH

#include "blade/base.hh"
#include "blade/module.hh"

namespace Blade::Modules {

template<typename IT, typename OT>
class BLADE_API Cast : public Module {
 public:
    // Configuration

    struct Config {
        U64 blockSize = 512;
    };

    constexpr const Config& getConfig() const {
        return this->config;
    }

    // Input

    struct Input {
        const ArrayTensor<Device::CUDA, IT>& buf;
    };

    constexpr const ArrayTensor<Device::CUDA, IT>& getInputBuffer() const {
        return this->input.buf;
    }

    // Output 

    struct Output {
        ArrayTensor<Device::CUDA, OT> buf;
    };

    constexpr const ArrayTensor<Device::CUDA, OT>& getOutputBuffer() const {
        return this->output.buf;
    }

    // Constructor & Processing

    explicit Cast(const Config& config, const Input& input);
    const Result process(const cudaStream_t& stream) final;

 private:
    // Variables

    const Config config;
    const Input input;
    Output output;

    // Expected Dimensions

    const ArrayDimensions getOutputBufferDims() const {
        return getInputBuffer().dims();
    }
};

}  // namespace Blade::Modules

#endif
