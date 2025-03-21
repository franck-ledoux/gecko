inc_dir = include_directories('inc')
lib_src = files('src/Blocking.cpp',
                'src/BlockingClassifier.cpp',
                'src/Graph.cpp')

gmdsblock_lib = library('gmdsblock_lib',
                    lib_src,
                    include_directories: inc_dir,
                    dependencies : [gmds_core_cad_lib_dep,
                                    gmds_core_cadfac_lib_dep,
                                    gmds_core_ig_lib_dep,
                                    gmds_core_igalgo_lib_dep,
                                    gmds_core_io_lib_dep,
                                    gmds_core_math_lib_dep,
                                    gmds_core_utils_lib_dep,
                                    eigen3_dep]
)

gmdsblock_lib_dep = declare_dependency(include_directories: inc_dir, link_with : gmdsblock_lib)

subdir('tst')