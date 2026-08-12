add_library(gmds_core_cadfac_lib
        src/FACCurve.cpp
        src/FACManager.cpp
        src/FACPoint.cpp
        src/FACSurface.cpp
        src/FACVolume.cpp
)
target_include_directories(gmds_core_cadfac_lib PUBLIC inc)
target_link_libraries(gmds_core_cadfac_lib PUBLIC
    gmds_core_cad_lib
    gmds_core_math_lib
    gmds_core_utils_lib
    gmds_core_ig_lib
    gmds_core_igalgo_lib
    gmds_core_io_lib
    boost_dep
    cgal_dep
    Eigen3::Eigen
)

add_subdirectory(tst)
