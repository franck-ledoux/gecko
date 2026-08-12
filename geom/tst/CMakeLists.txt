add_executable(test_geometry geometry_test.cpp)
target_link_libraries(test_geometry PRIVATE
    Gecko::Geom
    Catch2::Catch2
    Catch2::Catch2WithMain
)
# To get access to the generated file "unit_test_config.h"
target_include_directories(test_geometry PRIVATE ${CMAKE_BINARY_DIR})

add_test(NAME test_geometry COMMAND test_geometry)
