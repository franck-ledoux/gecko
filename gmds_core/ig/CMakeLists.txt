add_library(gmds_core_ig_lib
    src/Cell.cpp
    src/Edge.cpp
    src/Face.cpp
    src/Node.cpp
    src/Region.cpp
    src/EdgeContainer.cpp
    src/FaceContainer.cpp
    src/NodeContainer.cpp
    src/RegionContainer.cpp
    src/Mesh.cpp
    src/MeshDoctor.cpp
)
target_include_directories(gmds_core_ig_lib PUBLIC inc)
target_link_libraries(gmds_core_ig_lib PUBLIC
    Eigen3::Eigen
    gmds_core_utils_lib
    gmds_core_math_lib
)

add_subdirectory(tst)
