add_library(gmds_core_io_lib
    src/GMSHWriter.cpp
    src/IGMeshIOService.cpp
    src/IMeshIOService.cpp
    src/IReader.cpp
    src/IWriter.cpp
    src/MeditReader.cpp
    src/VTKReader.cpp
    src/VTKWriter.cpp
)
target_include_directories(gmds_core_io_lib PUBLIC inc)
target_link_libraries(gmds_core_io_lib PUBLIC
    Eigen3::Eigen
    gmds_core_utils_lib
    gmds_core_ig_lib
    gmds_core_math_lib
)

add_subdirectory(tst)
