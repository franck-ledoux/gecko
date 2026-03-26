add_executable(test_gmds_core_math
    math_tests.cpp
    point_tests.cpp
)
target_link_libraries(test_gmds_core_math PRIVATE
    gmds_core_math_lib
    gmds_core_utils_lib
    Catch2::Catch2
    Catch2::Catch2WithMain
)
add_test(NAME test_gmds_core_math COMMAND test_gmds_core_math)
