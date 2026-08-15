#include <gecko/math/ClosestPoint.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
using Catch::Approx;
using namespace gecko;

constexpr double EPSILON = 1e-9;

// =============================================================================
// closest_point_on_segment TESTS
// =============================================================================

TEST_CASE("closest_point_on_segment regions", "[math][ClosestPoint][Segment]") {
    const Point3d a(0.0, 0.0, 0.0);
    const Point3d b(10.0, 0.0, 0.0);

    SECTION("Point projects inside the segment") {
        const Point3d p(4.0, 3.0, 0.0);
        const Point3d cp = closest_point_on_segment(p, a, b);
        REQUIRE(cp.x() == Approx(4.0).margin(EPSILON));
        REQUIRE(cp.y() == Approx(0.0).margin(EPSILON));
    }

    SECTION("Point projects before the segment start (clamped to a)") {
        const Point3d p(-5.0, 2.0, 0.0);
        const Point3d cp = closest_point_on_segment(p, a, b);
        REQUIRE(cp == a);
    }

    SECTION("Point projects past the segment end (clamped to b)") {
        const Point3d p(15.0, -2.0, 0.0);
        const Point3d cp = closest_point_on_segment(p, a, b);
        REQUIRE(cp == b);
    }

    SECTION("Point already lies on the segment") {
        const Point3d p(6.0, 0.0, 0.0);
        const Point3d cp = closest_point_on_segment(p, a, b);
        REQUIRE(cp == p);
    }

    SECTION("Degenerate zero-length segment returns its single point") {
        const Point3d single(2.0, 2.0, 2.0);
        const Point3d p(9.0, -1.0, 4.0);
        const Point3d cp = closest_point_on_segment(p, single, single);
        REQUIRE(cp == single);
    }
}

// =============================================================================
// closest_point_on_triangle TESTS
// =============================================================================

TEST_CASE("closest_point_on_triangle Voronoi regions", "[math][ClosestPoint][Triangle]") {
    // Right triangle in the z = 0 plane.
    const Point3d a(0.0, 0.0, 0.0);
    const Point3d b(4.0, 0.0, 0.0);
    const Point3d c(0.0, 4.0, 0.0);

    SECTION("Point above the interior projects onto the face") {
        const Point3d p(1.0, 1.0, 5.0);
        const Point3d cp = closest_point_on_triangle(p, a, b, c);
        REQUIRE(cp.x() == Approx(1.0).margin(EPSILON));
        REQUIRE(cp.y() == Approx(1.0).margin(EPSILON));
        REQUIRE(cp.z() == Approx(0.0).margin(EPSILON));
    }

    SECTION("Point closest to vertex a") {
        const Point3d p(-3.0, -3.0, 0.0);
        const Point3d cp = closest_point_on_triangle(p, a, b, c);
        REQUIRE(cp == a);
    }

    SECTION("Point closest to vertex b") {
        const Point3d p(8.0, -3.0, 0.0);
        const Point3d cp = closest_point_on_triangle(p, a, b, c);
        REQUIRE(cp == b);
    }

    SECTION("Point closest to vertex c") {
        const Point3d p(-3.0, 8.0, 0.0);
        const Point3d cp = closest_point_on_triangle(p, a, b, c);
        REQUIRE(cp == c);
    }

    SECTION("Point closest to edge ab") {
        const Point3d p(2.0, -3.0, 0.0);
        const Point3d cp = closest_point_on_triangle(p, a, b, c);
        REQUIRE(cp.x() == Approx(2.0).margin(EPSILON));
        REQUIRE(cp.y() == Approx(0.0).margin(EPSILON));
    }

    SECTION("Point closest to edge ac") {
        const Point3d p(-3.0, 2.0, 0.0);
        const Point3d cp = closest_point_on_triangle(p, a, b, c);
        REQUIRE(cp.x() == Approx(0.0).margin(EPSILON));
        REQUIRE(cp.y() == Approx(2.0).margin(EPSILON));
    }

    SECTION("Point closest to edge bc") {
        const Point3d p(4.0, 4.0, 0.0);
        const Point3d cp = closest_point_on_triangle(p, a, b, c);
        // bc is the line x + y = 4; the closest point to (4,4) on it is (2,2).
        REQUIRE(cp.x() == Approx(2.0).margin(EPSILON));
        REQUIRE(cp.y() == Approx(2.0).margin(EPSILON));
    }

    SECTION("Point already on the triangle's plane and inside it") {
        const Point3d p(1.0, 1.0, 0.0);
        const Point3d cp = closest_point_on_triangle(p, a, b, c);
        REQUIRE(cp == p);
    }

    SECTION("Degenerate zero-area (collinear) triangle falls back to the closest edge") {
        const Point3d da(0.0, 0.0, 0.0);
        const Point3d db(2.0, 0.0, 0.0);
        const Point3d dc(4.0, 0.0, 0.0); // collinear with da, db
        const Point3d p(1.0, 3.0, 0.0);
        const Point3d cp = closest_point_on_triangle(p, da, db, dc);
        REQUIRE(cp.x() == Approx(1.0).margin(EPSILON));
        REQUIRE(cp.y() == Approx(0.0).margin(EPSILON));
    }
}
