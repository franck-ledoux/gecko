add_library(block_lib
    src/Blocking.cpp
    src/BlockingClassifier.cpp
    src/Graph.cpp
)
target_include_directories(block_lib PUBLIC inc)
target_link_libraries(block_lib PUBLIC
    cgal_dep
    boost_dep
    gmds_core_cad_lib
    gmds_core_cadfac_lib
    gmds_core_ig_lib
    gmds_core_igalgo_lib
    gmds_core_io_lib
    gmds_core_math_lib
    gmds_core_utils_lib
    Eigen3::Eigen
)

add_subdirectory(tst)
