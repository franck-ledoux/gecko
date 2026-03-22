add_executable(test_gmds_core_igalgo
    boundary_extractor_tests.cpp
    boundary_operator_tests.cpp
    grid_builder_tests.cpp
    surface_reorient_tests.cpp
)
target_link_libraries(test_gmds_core_igalgo PRIVATE
    Eigen3::Eigen
    gmds_core_math_lib
    gmds_core_io_lib
    gmds_core_ig_lib
    gmds_core_igalgo_lib
    gmds_core_utils_lib
    Catch2::Catch2
    Catch2::Catch2WithMain
    ut_dep
)
add_test(NAME test_gmds_core_igalgo COMMAND test_gmds_core_igalgo)
