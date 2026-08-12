add_executable(test_gmds_core_cadfac
        fac_manager_tests.cpp
        geom_linker_tests.cpp
        geom_topology_tests.cpp
)
target_link_libraries(test_gmds_core_cadfac PRIVATE
    gmds_core_cad_lib
    gmds_core_cadfac_lib
    gmds_core_math_lib
    gmds_core_utils_lib
    gmds_core_ig_lib
    gmds_core_igalgo_lib
    gmds_core_io_lib
    boost_dep
    cgal_dep
    Eigen3::Eigen
    Catch2::Catch2
    Catch2::Catch2WithMain
    ut_dep
)
add_test(NAME test_gmds_core_cadfac COMMAND test_gmds_core_cadfac)
