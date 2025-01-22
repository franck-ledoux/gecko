inc_dir = include_directories('inc')
lib_src = files('src/BitVector.cpp',
                'src/CommonTypes.cpp',
                'src/VariableManager.cpp')

gmds_core_utils_lib = library('gmds_core_utils_lib',
                     lib_src,
                     include_directories : inc_dir
)

gmds_core_utils_lib_dep = declare_dependency(include_directories : inc_dir, link_with : gmds_core_utils_lib)

subdir('tst')