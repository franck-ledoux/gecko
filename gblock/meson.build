inc_dir = include_directories('inc')
lib_src = files('src/Blocking.cpp')

gblock_lib = library('gblock_lib',
                    lib_src,
                    include_directories: inc_dir,
                     dependencies : [cgal_dep]
)

gblock_lib_dep = declare_dependency(include_directories: inc_dir, link_with : gblock_lib)