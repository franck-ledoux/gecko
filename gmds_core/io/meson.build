inc_dir = include_directories('inc')
lib_src = files('src/GMSHWriter.cpp',
                'src/IGMeshIOService.cpp',
                'src/IMeshIOService.cpp',
                'src/IReader.cpp',
                'src/IWriter.cpp',
                'src/MeditReader.cpp',
                'src/VTKReader.cpp',
                'src/VTKWriter.cpp')

gmds_core_io_lib = library('gmds_core_io_lib',
                     lib_src,
                     include_directories : inc_dir,
                     dependencies : [eigen3_dep, gmds_core_utils_lib_dep, gmds_core_ig_lib_dep, gmds_core_math_lib_dep]
)

gmds_core_io_lib_dep = declare_dependency(include_directories : inc_dir, link_with : gmds_core_io_lib)

subdir('tst')