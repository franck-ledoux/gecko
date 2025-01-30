project('predicates', 'cpp',
        version: '1.0.0')

inc_dir = include_directories('inc')
lib_src = files('src/Predicates_psm.cpp')

predicates_lib = library('ext_predicates_lib',
                     lib_src,
                     include_directories : inc_dir
)

predicates_lib_dep = declare_dependency(include_directories : inc_dir, link_with : predicates_lib)

