tst_src = files('reader_tests.cpp', 'writer_tests.cpp')

tst_exe = executable('test_gmds_core_math',
                      tst_src,
                      dependencies : [gmds_core_math_lib_dep,
                                      gmds_core_utils_lib_dep,
                                      gmds_core_ig_lib_dep,
                                      gmds_core_io_lib_dep,
                                      catch2_dep,
                                      catch2_main_dep ])

test('test_gmds_core_math', tst_exe)