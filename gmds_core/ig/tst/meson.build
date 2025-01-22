tst_src = files('boolean_mark_tests.cpp', 'doctor_tests.cpp','edge_tests.cpp','mesh_tests.cpp')
tst_exe = executable('test_gmds_core_ig',
                      tst_src,
                      dependencies : [eigen3_dep, gmds_core_math_lib_dep, gmds_core_ig_lib_dep, gmds_core_utils_lib_dep, catch2_dep, catch2_main_dep ])

test('test_gmds_core_ig', tst_exe)