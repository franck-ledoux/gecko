add_executable(test_utils utils_tests.cpp)
target_link_libraries(test_utils PRIVATE
    Gecko::Utils
    Catch2::Catch2
    Catch2::Catch2WithMain
)
# To get access to the generated file "unit_test_config.h"
target_include_directories(test_utils PRIVATE ${CMAKE_BINARY_DIR})

add_test(NAME test_utils COMMAND test_utils)
