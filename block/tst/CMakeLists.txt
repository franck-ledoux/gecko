add_executable(test_cblock cblock_tests.cpp)
target_link_libraries(test_cblock PRIVATE
    Eigen3::Eigen
    cgal_dep
    boost_dep
    gmds_core_cad_lib
    gmds_core_cadfac_lib
    gmds_core_ig_lib
    gmds_core_igalgo_lib
    gmds_core_io_lib
    gmds_core_math_lib
    gmds_core_utils_lib
    cblock_lib
    Catch2::Catch2
    Catch2::Catch2WithMain

)
# To get access to the generated file "unit_test_config.h"
target_include_directories(test_cblock PRIVATE ${CMAKE_BINARY_DIR})

add_test(NAME test_cblock COMMAND test_cblock)
