add_library(gmds_core_cad_lib
        src/GeomSurface.cpp
        src/GeomMeshLinker.cpp
)
target_include_directories(gmds_core_cad_lib PUBLIC inc)
target_link_libraries(gmds_core_cad_lib PUBLIC
    gmds_core_math_lib
    gmds_core_io_lib
    gmds_core_utils_lib
    gmds_core_ig_lib
    Eigen3::Eigen
)
