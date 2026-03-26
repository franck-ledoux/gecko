add_executable(test_gmds_core_utils utils_tests.cpp)
target_link_libraries(test_gmds_core_utils PRIVATE
    gmds_core_utils_lib
    Catch2::Catch2
    Catch2::Catch2WithMain
)
add_test(NAME test_gmds_core_utils COMMAND test_gmds_core_utils)
