deps = [
]

all_deps_found = true
foreach x_dep : deps
    all_deps_found = all_deps_found and cfg_lst.get(x_dep, false)
endforeach

if all_deps_found
    src_lst += files([
        'base.cc',
    ])

    jit_lst += files([
        'detector.cu',
    ])
endif

sum_lst += {'Detector': all_deps_found}

cfg_lst.set('BLADE_MODULE_DETECTOR', all_deps_found)
