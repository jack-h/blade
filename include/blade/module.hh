#ifndef BLADE_MODULE_HH
#define BLADE_MODULE_HH

#include <map>
#include <string>
#include <typeindex>

#include "blade/types.hh"
#include "blade/logger.hh"
#include "blade/macros.hh"

#include "blade/utils/jitify2.hh"
using namespace jitify2::reflection;

namespace Blade {

class Module {
 public:
    explicit Module(const U64& blockSize,
                    const jitify2::PreprocessedProgram& kernel);
    virtual ~Module() = default;

    virtual constexpr const Result  preprocess(const cudaStream_t& stream = 0) {
        return Result::SUCCESS;
    }

    virtual constexpr const Result process(const cudaStream_t& stream = 0) {
        return Result::SUCCESS;
    }

 protected:
    jitify2::ProgramCache<> cache;
    std::string kernel;
    dim3 grid, block;

    template<typename T>
    static const std::string CudaType() {
        static std::unordered_map<std::type_index, std::string> type_map = {
            {typeid(CF16),  "__half"},
            {typeid(CF32),  "float"},
            {typeid(CI8),   "char"},
            {typeid(CI16),  "short"},
            {typeid(CI32),  "long"},
            {typeid(F16),   "__half"},
            {typeid(F32),   "float"},
            {typeid(I8),    "char"},
            {typeid(I16),   "short"},
            {typeid(I32),   "long"},
        };

        auto& type = typeid(T);
        if (!type_map.contains(type)) {
            BL_FATAL("Type is not supported by CudaType.");
            BL_CHECK_THROW(Result::ERROR);
        }
        return type_map[type];
    }

    template<typename T>
    static const std::size_t CudaTypeSize() {
        static std::unordered_map<std::type_index, std::size_t> size_map = {
            {typeid(CF16),  2},
            {typeid(CF32),  2},
            {typeid(CI8),   2},
            {typeid(CI16),  2},
            {typeid(CI32),  2},
            {typeid(F16),   1},
            {typeid(F32),   1},
            {typeid(I8),    1},
            {typeid(I16),   1},
            {typeid(I32),   1},
        };

        auto& type = typeid(T);
        if (!size_map.contains(type)) {
            BL_FATAL("Type is not supported by CudaTypeSize.");
            BL_CHECK_THROW(Result::ERROR);
        }
        return size_map[type];
    }
};

}  // namespace Blade

#endif
