#include <gecko/math/BezierCurve.h>
#include <gecko/math/Point3d.h>
#include <gecko/math/Vector3d.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
using Catch::Approx;
using namespace gecko;

TEST_CASE("bezier_curve_degree_1_is_a_straight_line", "[MathTestSuite]") {
    const BezierCurve<1, Point3d> line(Point3d(0.0, 0.0, 0.0), Point3d(4.0, 0.0, 0.0));

    REQUIRE(line.value(0.0) == Point3d(0.0, 0.0, 0.0));
    REQUIRE(line.value(1.0) == Point3d(4.0, 0.0, 0.0));
    const Point3d mid = line.value(0.5);
    REQUIRE(mid.x() == Approx(2.0));
    REQUIRE(mid.y() == Approx(0.0));
    REQUIRE(mid.z() == Approx(0.0));
}

TEST_CASE("bezier_curve_degree_2_evaluates_via_de_casteljau", "[MathTestSuite]") {
    // Quadratic Bezier with a control point pulled off the chord: bulges toward it.
    const BezierCurve<2, Point3d> curve(Point3d(0.0, 0.0, 0.0), Point3d(2.0, 4.0, 0.0), Point3d(4.0, 0.0, 0.0));

    REQUIRE(curve.value(0.0) == Point3d(0.0, 0.0, 0.0));
    REQUIRE(curve.value(1.0) == Point3d(4.0, 0.0, 0.0));
    const Point3d mid = curve.value(0.5);
    // B(0.5) = 0.25*P0 + 0.5*P1 + 0.25*P2 = (2.0, 2.0, 0.0)
    REQUIRE(mid.x() == Approx(2.0));
    REQUIRE(mid.y() == Approx(2.0));
    REQUIRE(mid.z() == Approx(0.0));
}

TEST_CASE("bezier_curve_default_constructed_control_points_are_origin", "[MathTestSuite]") {
    const BezierCurve<2, Point3d> curve;
    REQUIRE(curve[0] == Point3d(0.0, 0.0, 0.0));
    REQUIRE(curve[1] == Point3d(0.0, 0.0, 0.0));
    REQUIRE(curve[2] == Point3d(0.0, 0.0, 0.0));
}
