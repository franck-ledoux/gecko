tst_src = files('utils_tests.cpp')

tst_exe = executable('test_gmds_core_utils',
                      tst_src,
                      dependencies : [ gmds_core_utils_lib_dep, catch2_dep, catch2_main_dep ])

test('test_gmds_core_utils', tst_exe)