add_executable(test_gmds_core_ig
    boolean_mark_tests.cpp
    doctor_tests.cpp
    edge_tests.cpp
    mesh_tests.cpp
)
target_link_libraries(test_gmds_core_ig PRIVATE
    Eigen3::Eigen
    gmds_core_math_lib
    gmds_core_ig_lib
    gmds_core_utils_lib
    Catch2::Catch2
    Catch2::Catch2WithMain
)
add_test(NAME test_gmds_core_ig COMMAND test_gmds_core_ig)
