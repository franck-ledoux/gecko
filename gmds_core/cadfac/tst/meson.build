tst_src = files('fac_manager_tests.cpp', 'geom_linker_tests.cpp', 'geom_topology_tests.cpp')
tst_exe = executable('test_gmds_core_cadfac',
                      tst_src,
                      dependencies : [gmds_core_cad_lib_dep,
                                      gmds_core_cadfac_lib_dep,
                                      gmds_core_math_lib_dep,
                                      gmds_core_utils_lib_dep,
                                      gmds_core_ig_lib_dep,
                                      gmds_core_igalgo_lib_dep,
                                      gmds_core_io_lib_dep,
                                      boost_dep,
                                      cgal_dep,
                                      eigen3_dep,
                                      catch2_dep,
                                      catch2_main_dep,
                                      ut_dep ])

test('test_gmds_core_cadfac', tst_exe)