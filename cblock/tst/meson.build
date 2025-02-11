tst_src = files('cblock_tests.cpp')

tst_exe = executable('test_cblock',
                      tst_src,
                      dependencies : [eigen3_dep,
                                      cgal_dep,
                                      boost_dep,
                                      gmds_core_cad_lib_dep,
                                      gmds_core_cadfac_lib_dep,
                                      gmds_core_ig_lib_dep,
                                      gmds_core_igalgo_lib_dep,
                                      gmds_core_io_lib_dep,
                                      gmds_core_math_lib_dep,
                                      gmds_core_utils_lib_dep,
                                      cblock_lib_dep,
                                      catch2_dep,
                                      catch2_main_dep,
                                      ut_dep ])

test('test_cblock', tst_exe)