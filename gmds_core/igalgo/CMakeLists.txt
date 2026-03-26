add_library(gmds_core_igalgo_lib
    src/BoundaryExtractor2D.cpp
    src/BoundaryExtractor3D.cpp
    src/BoundaryOperator.cpp
    src/BoundaryOperator2D.cpp
    src/GridBuilder.cpp
    src/SurfaceReorient.cpp
    src/THexBuilder.cpp
)
target_include_directories(gmds_core_igalgo_lib PUBLIC inc)
target_link_libraries(gmds_core_igalgo_lib PUBLIC
    Eigen3::Eigen
    gmds_core_ig_lib
    gmds_core_utils_lib
    gmds_core_math_lib
)

add_subdirectory(tst)
