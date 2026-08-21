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

TEST_CASE("bezier_curve_split_reproduces_the_parent_exactly", "[MathTestSuite]") {
    // The property the whole block-structure cut rests on: subdivision is exact, not a fit. A
    // deliberately asymmetric cubic, so a wrong-but-plausible implementation (mirrored halves,
    // swapped pyramid edges) cannot pass by symmetry.
    const BezierCurve<3, Point3d> curve(
        Point3d(0.0, 0.0, 0.0), Point3d(1.0, 3.0, 0.5), Point3d(4.0, -1.0, 2.0), Point3d(5.0, 2.0, 1.0));

    const double s = 0.37;
    const auto [left, right] = curve.split(s);

    for (int i = 0; i <= 20; ++i) {
        const double t = static_cast<double>(i) / 20.0;
        const Point3d on_left = left.value(t);
        const Point3d expected_left = curve.value(t * s);
        REQUIRE(on_left.x() == Approx(expected_left.x()).margin(1e-12));
        REQUIRE(on_left.y() == Approx(expected_left.y()).margin(1e-12));
        REQUIRE(on_left.z() == Approx(expected_left.z()).margin(1e-12));

        const Point3d on_right = right.value(t);
        const Point3d expected_right = curve.value(s + t * (1.0 - s));
        REQUIRE(on_right.x() == Approx(expected_right.x()).margin(1e-12));
        REQUIRE(on_right.y() == Approx(expected_right.y()).margin(1e-12));
        REQUIRE(on_right.z() == Approx(expected_right.z()).margin(1e-12));
    }

    // The 2 halves must meet exactly at the split point, so a node placed there is unambiguous.
    REQUIRE(left.control_points()[3] == right.control_points()[0]);
    REQUIRE(left.control_points()[0] == curve.control_points()[0]);
    REQUIRE(right.control_points()[3] == curve.control_points()[3]);
}

TEST_CASE("bezier_curve_split_of_a_straight_line_stays_straight", "[MathTestSuite]") {
    // Degree 1 is the linear-blocking case, and its halves must stay exact segments of the parent.
    const BezierCurve<1, Point3d> line(Point3d(0.0, 0.0, 0.0), Point3d(4.0, 8.0, 0.0));
    const auto [left, right] = line.split(0.25);

    REQUIRE(left.control_points()[0] == Point3d(0.0, 0.0, 0.0));
    REQUIRE(left.control_points()[1].x() == Approx(1.0));
    REQUIRE(left.control_points()[1].y() == Approx(2.0));
    REQUIRE(right.control_points()[0] == left.control_points()[1]);
    REQUIRE(right.control_points()[1] == Point3d(4.0, 8.0, 0.0));
}
