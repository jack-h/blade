#include "blade/beamformer/test/generic.hh"
#include "blade/beamformer/generic.hh"
#include "blade/checker/base.hh"
#include "blade/manager.hh"
#include "blade/pipeline.hh"

using namespace Blade;

template<typename T, typename TT>
class Module : public Pipeline {
public:
    Module(const typename T::Config &config) : config(config) {
        if (this->commit() != Result::SUCCESS) {
            throw Result::ERROR;
        }
    }

protected:
    Result underlyingInit() final {
        BL_INFO("Initializing kernels.");

        beamformer = Factory<T>(config);
        test = std::make_unique<TT>(config);
        checker = Factory<Checker::Generic>({});

        return Result::SUCCESS;
    }

    Result underlyingAllocate() final {
        BL_INFO("Allocating resources.");

        BL_CHECK(allocateBuffer(input, beamformer->getInputSize(), true));
        BL_CHECK(allocateBuffer(phasors, beamformer->getPhasorsSize(), true));
        BL_CHECK(allocateBuffer(output, beamformer->getOutputSize(), true));
        BL_CHECK(allocateBuffer(result, beamformer->getOutputSize(), true));

        BL_INFO("Generating test data.");
        for (auto & element : input) {
            element = 1;
        }

        for (auto & element : phasors) {
            element = 2.1;
        }

        for (auto & element : result) {
            element = config.NANTS * 2.1;
        }

        BL_CHECK(test->process());

        return Result::SUCCESS;
    }

    Result underlyingReport(Resources &res) final {
        BL_INFO("Reporting resources.");

        res.transfer.h2d += input.size_bytes();
        res.transfer.h2d += phasors.size_bytes();
        res.transfer.d2h += output.size_bytes();

        return Result::SUCCESS;
    }

    Result underlyingProcess(cudaStream_t &cudaStream) final {
        BL_CHECK(beamformer->run(input, phasors, output, cudaStream));

        return Result::SUCCESS;
    }

    Result underlyingPostprocess() final {
        std::size_t errors = 0;
        if ((errors = checker->run(output, result)) != 0) {
            BL_FATAL("Module produced {} errors.", errors);
            return Result::ERROR;
        }

        return Result::SUCCESS;
    }

private:
    const typename T::Config &config;

    std::unique_ptr<Beamformer::Generic> beamformer;
    std::unique_ptr<Beamformer::Test::Generic> test;
    std::unique_ptr<Checker::Generic> checker;

    std::span<std::complex<float>> input;
    std::span<std::complex<float>> phasors;
    std::span<std::complex<float>> output;
    std::span<std::complex<float>> result;
};
