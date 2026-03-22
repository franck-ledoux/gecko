add_executable(test_gmds_core_io
    reader_tests.cpp
    writer_tests.cpp
)
target_link_libraries(test_gmds_core_io PRIVATE
    gmds_core_math_lib
    gmds_core_utils_lib
    gmds_core_ig_lib
    gmds_core_io_lib
    Eigen3::Eigen
    Catch2::Catch2
    Catch2::Catch2WithMain
    ut_dep
)
add_test(NAME test_gmds_core_io COMMAND test_gmds_core_io)
