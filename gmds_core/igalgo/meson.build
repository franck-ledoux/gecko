inc_dir = include_directories('inc')
lib_src = files('src/BoundaryExtractor2D.cpp',
                'src/BoundaryExtractor3D.cpp',
                'src/BoundaryOperator.cpp',
                'src/BoundaryOperator2D.cpp',
                'src/GridBuilder.cpp',
                'src/SurfaceReorient.cpp',
                'src/THexBuilder.cpp')

gmds_core_igalgo_lib = library('gmds_core_igalgo_lib',
                     lib_src,
                     include_directories : inc_dir,
                     dependencies : [eigen3_dep, gmds_core_ig_lib_dep, gmds_core_utils_lib_dep, gmds_core_math_lib_dep]
)

gmds_core_igalgo_lib_dep = declare_dependency(include_directories : inc_dir, link_with : gmds_core_igalgo_lib)

subdir('tst')