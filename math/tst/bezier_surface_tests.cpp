#include <gecko/math/BezierSurface.h>
#include <gecko/math/Point3d.h>
#include <gecko/math/Vector3d.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
using Catch::Approx;
using namespace gecko;

TEST_CASE("bezier_surface_degree_1_is_bilinear", "[MathTestSuite]") {
    // Flat unit square in the z=0 plane: corners [i][j], i along u, j along v.
    BezierSurface<1, Point3d>::ControlGrid grid;
    grid[0][0] = Point3d(0.0, 0.0, 0.0); // u=0,v=0
    grid[0][1] = Point3d(0.0, 1.0, 0.0); // u=0,v=1
    grid[1][0] = Point3d(1.0, 0.0, 0.0); // u=1,v=0
    grid[1][1] = Point3d(1.0, 1.0, 0.0); // u=1,v=1
    const BezierSurface<1, Point3d> surface(grid);

    REQUIRE(surface.value(0.0, 0.0) == Point3d(0.0, 0.0, 0.0));
    REQUIRE(surface.value(1.0, 0.0) == Point3d(1.0, 0.0, 0.0));
    REQUIRE(surface.value(0.0, 1.0) == Point3d(0.0, 1.0, 0.0));
    REQUIRE(surface.value(1.0, 1.0) == Point3d(1.0, 1.0, 0.0));

    const Point3d center = surface.value(0.5, 0.5);
    REQUIRE(center.x() == Approx(0.5));
    REQUIRE(center.y() == Approx(0.5));
    REQUIRE(center.z() == Approx(0.0));
}

TEST_CASE("bezier_surface_degree_2_bulges_toward_interior_control_point", "[MathTestSuite]") {
    BezierSurface<2, Point3d>::ControlGrid grid;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            grid[i][j] = Point3d(static_cast<double>(i), static_cast<double>(j), 0.0);
        }
    }
    // Pull the center control point up out of the plane.
    grid[1][1] = Point3d(1.0, 1.0, 4.0);
    const BezierSurface<2, Point3d> surface(grid);

    // Corners are still interpolated exactly.
    REQUIRE(surface.value(0.0, 0.0) == Point3d(0.0, 0.0, 0.0));
    REQUIRE(surface.value(1.0, 1.0) == Point3d(2.0, 2.0, 0.0));

    // The surface center bulges toward the raised control point (z > 0).
    const Point3d center = surface.value(0.5, 0.5);
    REQUIRE(center.z() > 0.5);
}

TEST_CASE("bezier_surface_default_constructed_control_points_are_origin", "[MathTestSuite]") {
    const BezierSurface<1, Point3d> surface;
    REQUIRE(surface.control_point(0, 0) == Point3d(0.0, 0.0, 0.0));
    REQUIRE(surface.control_point(1, 1) == Point3d(0.0, 0.0, 0.0));
}
