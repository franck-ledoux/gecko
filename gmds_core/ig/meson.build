inc_dir = include_directories('inc')
lib_src = files('src/Cell.cpp',
                'src/Edge.cpp',
                'src/Face.cpp',
                'src/Node.cpp',
                'src/Region.cpp',
                'src/EdgeContainer.cpp',
                'src/FaceContainer.cpp',
                'src/NodeContainer.cpp',
                'src/RegionContainer.cpp',
                'src/Mesh.cpp',
                'src/MeshDoctor.cpp')

gmds_core_ig_lib = library('gmds_core_ig_lib',
                     lib_src,
                     include_directories : inc_dir,
                     dependencies : [eigen3_dep, gmds_core_utils_lib_dep, gmds_core_math_lib_dep]
)

gmds_core_ig_lib_dep = declare_dependency(include_directories : inc_dir, link_with : gmds_core_ig_lib)

subdir('tst')