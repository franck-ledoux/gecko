add_executable(test_mesh mesh_tests.cpp)
target_link_libraries(test_mesh PRIVATE
    Gecko::Mesh
    Catch2::Catch2
    Catch2::Catch2WithMain
)
# To get access to the generated file "unit_test_config.h"
target_include_directories(test_mesh PRIVATE ${CMAKE_BINARY_DIR})

add_test(NAME test_mesh COMMAND test_mesh)
