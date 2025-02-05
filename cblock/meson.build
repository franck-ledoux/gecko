inc_dir = include_directories('inc')
lib_src = files('src/Blocking.cpp')

cblock_lib = library('cblock_lib',
                    lib_src,
                    include_directories: inc_dir,
                    dependencies : [cgal_dep,
                                    boost_dep,
                                    gmds_core_cad_lib_dep,
                                    gmds_core_cadfac_lib_dep,
                                    gmds_core_ig_lib_dep,
                                    gmds_core_igalgo_lib_dep,
                                    gmds_core_io_lib_dep,
                                    gmds_core_math_lib_dep,
                                    gmds_core_utils_lib_dep,
                                    eigen3_dep]
)

cblock_lib_dep = declare_dependency(include_directories: inc_dir, link_with : cblock_lib)

subdir('tst')