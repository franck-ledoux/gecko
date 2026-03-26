add_library(gmds_core_math_lib
    src/Numerics.cpp
    src/AxisAngleRotation.cpp
    src/Chart.cpp
    src/Cross.cpp
    src/Cross2D.cpp
    src/Point.cpp
    src/Line.cpp
    src/Plane.cpp
    src/Orientation.cpp
    src/Quadrilateral.cpp
    src/Quaternion.cpp
    src/Segment.cpp
    src/Vector.cpp
    src/VectorDyn.cpp
    src/Hexahedron.cpp
    src/Prism3.cpp
    src/Pyramid.cpp
    src/Tetrahedron.cpp
    src/Triangle.cpp
    src/Ray.cpp
)
target_include_directories(gmds_core_math_lib PUBLIC inc)
target_link_libraries(gmds_core_math_lib PUBLIC
    Eigen3::Eigen
    ext_predicates_lib
    gmds_core_utils_lib
)

add_subdirectory(tst)
