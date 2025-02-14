inc_dir = include_directories('inc', 'spam/inc')

lib_src = files('src/BlockingAction.cpp',
                'src/BlockingState.cpp',
                'src/BlockingRewardFunction.cpp',
                'spam/src/MCTSAgent.cpp',
                'spam/src/MCTSLoop.cpp',
                'spam/src/MCTSSelectionFunction.cpp',
                'spam/src/MCTSTree.cpp')

exe_src = files('src/main_blocker.cpp')

mcts_c_lib = library('mcts_lib',
                    lib_src,
                    include_directories: inc_dir,
                    dependencies : [cgal_dep,
                                    boost_dep,
                                    nlohmann_json_dep,
                                    gmds_core_cad_lib_dep,
                                    gmds_core_cadfac_lib_dep,
                                    gmds_core_ig_lib_dep,
                                    gmds_core_igalgo_lib_dep,
                                    gmds_core_io_lib_dep,
                                    gmds_core_math_lib_dep,
                                    gmds_core_utils_lib_dep,
                                    cblock_lib_dep,
                                    eigen3_dep]
)

mcts_c_lib_dep = declare_dependency(include_directories: inc_dir, link_with : mcts_c_lib)

mcts_exe = executable('mcts_blocker', exe_src, dependencies: [mcts_c_lib_dep,cgal_dep,
                                                              boost_dep,
                                                              nlohmann_json_dep,
                                                              gmds_core_cad_lib_dep,
                                                              gmds_core_cadfac_lib_dep,
                                                              gmds_core_ig_lib_dep,
                                                              gmds_core_igalgo_lib_dep,
                                                              gmds_core_io_lib_dep,
                                                              gmds_core_math_lib_dep,
                                                              gmds_core_utils_lib_dep,
                                                              cblock_lib_dep,
                                                              eigen3_dep])