inc_dir = include_directories('inc')
lib_src = files('src/Numerics.cpp',
                'src/AxisAngleRotation.cpp',
                'src/Chart.cpp',
                'src/Cross.cpp',
                'src/Cross2D.cpp',
                'src/Point.cpp',
                'src/Line.cpp',
                'src/Plane.cpp',
                'src/Orientation.cpp',
                'src/Quadrilateral.cpp',
                'src/Quaternion.cpp',
                'src/Segment.cpp',
                'src/Vector.cpp',
                'src/VectorDyn.cpp',
                'src/Hexahedron.cpp',
                'src/Prism3.cpp',
                'src/Pyramid.cpp',
                'src/Tetrahedron.cpp',
                'src/Triangle.cpp',
                'src/Ray.cpp')

gmds_core_math_lib = library('gmds_core_math_lib',
                     lib_src,
                     include_directories : inc_dir,
                     dependencies : [eigen3_dep, ext_predicates_lib_dep, gmds_core_utils_lib_dep]
)

gmds_core_math_lib_dep = declare_dependency(include_directories : inc_dir, link_with : gmds_core_math_lib)

subdir('tst')