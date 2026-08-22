#include <gecko/math/BezierSurface.h>
#include <gecko/math/Point3d.h>
#include <gecko/math/Vector3d.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
using Catch::Approx;
using namespace gecko;

TEST_CASE("bezier_surface_degree_1_is_bilinear", "[MathTestSuite]") {
    // Flat unit square in the z=0 plane: corners [i][j], i along u, j along v.
    auto grid = BezierSurface<Point3d>::make_grid(1);
    grid[0][0] = Point3d(0.0, 0.0, 0.0); // u=0,v=0
    grid[0][1] = Point3d(0.0, 1.0, 0.0); // u=0,v=1
    grid[1][0] = Point3d(1.0, 0.0, 0.0); // u=1,v=0
    grid[1][1] = Point3d(1.0, 1.0, 0.0); // u=1,v=1
    const BezierSurface<Point3d> surface(grid);

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
    auto grid = BezierSurface<Point3d>::make_grid(2);
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            grid[i][j] = Point3d(static_cast<double>(i), static_cast<double>(j), 0.0);
        }
    }
    // Pull the center control point up out of the plane.
    grid[1][1] = Point3d(1.0, 1.0, 4.0);
    const BezierSurface<Point3d> surface(grid);

    // Corners are still interpolated exactly.
    REQUIRE(surface.value(0.0, 0.0) == Point3d(0.0, 0.0, 0.0));
    REQUIRE(surface.value(1.0, 1.0) == Point3d(2.0, 2.0, 0.0));

    // The surface center bulges toward the raised control point (z > 0).
    const Point3d center = surface.value(0.5, 0.5);
    REQUIRE(center.z() > 0.5);
}

TEST_CASE("bezier_surface_default_constructed_control_points_are_origin", "[MathTestSuite]") {
    const BezierSurface<Point3d> surface;
    REQUIRE(surface.control_point(0, 0) == Point3d(0.0, 0.0, 0.0));
    REQUIRE(surface.control_point(1, 1) == Point3d(0.0, 0.0, 0.0));
}

TEST_CASE("bezier_surface_split_u_and_split_v_reproduce_the_parent_exactly", "[MathTestSuite]") {
    // Same exactness property as BezierCurve::split(), one dimension up: a wrong axis (u/v swapped)
    // or a wrong fiber gather shows up immediately, since the control grid below is not symmetric.
    auto grid = BezierSurface<Point3d>::make_grid(2);
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            const double x = static_cast<double>(i);
            const double y = static_cast<double>(j) * 1.5;
            const double z = 0.4 * static_cast<double>(i * i) - 0.3 * static_cast<double>(j) + 0.2;
            grid[i][j] = Point3d(x, y, z);
        }
    }
    const BezierSurface<Point3d> surface(grid);

    const double su = 0.42;
    const auto [u_low, u_high] = surface.split_u(su);
    const double sv = 0.63;
    const auto [v_low, v_high] = surface.split_v(sv);

    for (int a = 0; a <= 6; ++a) {
        for (int b = 0; b <= 6; ++b) {
            const double u = static_cast<double>(a) / 6.0;
            const double v = static_cast<double>(b) / 6.0;

            const Point3d lu = u_low.value(u, v);
            const Point3d elu = surface.value(u * su, v);
            REQUIRE(lu.x() == Approx(elu.x()).margin(1e-12));
            REQUIRE(lu.y() == Approx(elu.y()).margin(1e-12));
            REQUIRE(lu.z() == Approx(elu.z()).margin(1e-12));

            const Point3d hu = u_high.value(u, v);
            const Point3d ehu = surface.value(su + u * (1.0 - su), v);
            REQUIRE(hu.z() == Approx(ehu.z()).margin(1e-12));
            REQUIRE(hu.x() == Approx(ehu.x()).margin(1e-12));

            const Point3d lv = v_low.value(u, v);
            const Point3d elv = surface.value(u, v * sv);
            REQUIRE(lv.y() == Approx(elv.y()).margin(1e-12));
            REQUIRE(lv.z() == Approx(elv.z()).margin(1e-12));

            const Point3d hv = v_high.value(u, v);
            const Point3d ehv = surface.value(u, sv + v * (1.0 - sv));
            REQUIRE(hv.y() == Approx(ehv.y()).margin(1e-12));
            REQUIRE(hv.z() == Approx(ehv.z()).margin(1e-12));
        }
    }
}
