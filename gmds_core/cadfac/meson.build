inc_dir = include_directories('inc')
lib_src = files('src/FACCurve.cpp',
                'src/FACManager.cpp',
                'src/FACPoint.cpp',
                'src/FACSurface.cpp',
                'src/FACVolume.cpp')

gmds_core_cadfac_lib = library('gmds_core_cadfac_lib',
                     lib_src,
                     include_directories : inc_dir,
                     dependencies : [gmds_core_cad_lib_dep,
                                     gmds_core_math_lib_dep,
                                     gmds_core_utils_lib_dep,
                                     gmds_core_ig_lib_dep,
                                     gmds_core_igalgo_lib_dep,
                                     gmds_core_io_lib_dep,
                                     eigen3_dep]
)

gmds_core_cadfac_lib_dep = declare_dependency(include_directories : inc_dir, link_with : gmds_core_cadfac_lib)

