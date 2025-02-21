project(
    'blade',
    ['cpp', 'c', 'cuda'],
    version: '0.7.1',
    default_options: [
        'buildtype=release',
        'cpp_std=c++20',
        'b_ndebug=if-release',
    ]
)

cfg_lst = configuration_data()
src_lst = []
dep_lst = [
    dependency('cuda', version : '>=11', modules : ['nvrtc', 'cudart', 'cuda', 'cufft']),
    dependency('fmt'),
]
inc_lst = [
    include_directories('include'),
    include_directories('.'),
]
jit_lst = []

subdir('meson')
subdir('src')
subdir('include')
subdir('tools')

lib_blade = library(
    'blade',
    src_lst,
    include_directories: inc_lst,
    dependencies: dep_lst,
    gnu_symbol_visibility: 'hidden',
    install: true,
)

lib_blade_dep = declare_dependency(
    include_directories: inc_lst,
    dependencies: dep_lst,
    link_with: lib_blade,
)

subdir('python')
subdir('tests')
subdir('docs')

summary({'cpp_std': get_option('cpp_std'),
            'prefix': get_option('prefix'),
            'buildtype': get_option('buildtype'),
            'python bindings': has_python,
        }, section: 'General', bool_yn: true)

pkg = import('pkgconfig')
pkg.generate(
    lib_blade,
    libraries: dep_lst,
)

subdir('apps')
