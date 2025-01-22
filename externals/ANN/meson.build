inc_dir = include_directories('include')
lib_src = files(
    'src/ANN.cpp',
    'src/bd_fix_rad_search.cpp',
    'src/bd_pr_search.cpp',
    'src/bd_search.cpp',
    'src/bd_tree.cpp',
    'src/bd_tree.h',
    'src/brute.cpp',
    'src/kd_dump.cpp',
    'src/kd_fix_rad_search.cpp',
    'src/kd_fix_rad_search.h',
    'src/kd_pr_search.cpp',
    'src/kd_pr_search.h',
    'src/kd_search.cpp',
    'src/kd_search.h',
    'src/kd_split.cpp',
    'src/kd_split.h',
    'src/kd_tree.cpp',
    'src/kd_tree.h',
    'src/kd_util.cpp',
    'src/kd_util.h',
    'src/perf.cpp',
    'src/pr_queue.h',
    'src/pr_queue_k.h')

ext_ann_lib = library('ext_ann_lib',
                     lib_src,
                     include_directories : inc_dir
)

ext_ann_lib_dep = declare_dependency(include_directories : inc_dir, link_with : ext_ann_lib)

