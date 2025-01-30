subdir('utils')
subdir('math')
subdir('ig')
subdir('io')
subdir('igalgo')
subdir('cad')
subdir('cadfac')

gmds_core_lib = library(
    'gmds_core_lib',
    [],
    dependencies: [gmds_core_cad_lib_dep,
                   gmds_core_cadfac_lib_dep,
                   gmds_core_ig_lib_dep,
                   gmds_core_igalgo_lib_dep,
                   gmds_core_io_lib_dep,
                   gmds_core_math_lib_dep,
                   gmds_core_utils_lib_dep,
                   cgal_dep,
                   eigen3_dep]

)


gmds_core_lib_dep = declare_dependency(include_directories : inc_dir, link_with : gmds_core_lib)
