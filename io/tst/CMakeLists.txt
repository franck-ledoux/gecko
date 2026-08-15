add_executable(test_io gmsh_io_tests.cpp)
target_link_libraries(test_io PRIVATE
    Gecko::Io
    Catch2::Catch2
    Catch2::Catch2WithMain
)
# To get access to the generated file "unit_test_config.h"
target_include_directories(test_io PRIVATE ${CMAKE_BINARY_DIR})

add_test(NAME test_io COMMAND test_io)
