inc_dir = include_directories('inc')
lib_src = files('src/Predicates_psm.cpp')

ext_predicates_lib = library('ext_predicates_lib',
                     lib_src,
                     include_directories : inc_dir
)

ext_predicates_lib_dep = declare_dependency(include_directories : inc_dir, link_with : ext_predicates_lib)

