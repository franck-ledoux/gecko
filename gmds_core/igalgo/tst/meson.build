tst_src = files('boundary_extractor_tests.cpp',
                'boundary_operator_tests.cpp',
                'grid_builder_tests.cpp',
                'surface_reorient_tests.cpp')

tst_exe = executable('test_gmds_core_igalgo',
                      tst_src,
                      dependencies : [eigen3_dep,
                                      gmds_core_math_lib_dep,
                                      gmds_core_io_lib_dep,
                                      gmds_core_ig_lib_dep,
                                      gmds_core_igalgo_lib_dep,
                                      gmds_core_utils_lib_dep,
                                      catch2_dep,
                                      catch2_main_dep,
                                      ut_dep ])

test('test_gmds_core_igalgo', tst_exe)