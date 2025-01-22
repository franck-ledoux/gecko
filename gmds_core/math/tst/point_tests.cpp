#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <sstream>
#include <string>
#include <gmds/math/Point.h>


using Catch::Approx;

TEST_CASE("PointClass - Setter") {
    gmds::math::Point p(1, 2, 3);
    gmds::TCoord z = p.Z();
    REQUIRE(p.X() == Approx(1).epsilon(1e-6));
    REQUIRE(p.Y() == Approx(2).epsilon(1e-6));
    REQUIRE(z == Approx(3).epsilon(1e-6));
}

TEST_CASE("PointClass - Setter2") {
    gmds::math::Point p(1, 2, 3);
    p.X() = 2 * p.X();
    p.Y() = 3 * p.Y();
    p.Z() = 0 * p.Z();

    REQUIRE(p.X() == Approx(2).epsilon(1e-6));
    REQUIRE(p.Y() == Approx(6).epsilon(1e-6));
    REQUIRE(p.Z() == Approx(0).epsilon(1e-6));
}

TEST_CASE("PointClass - SetXYZ") {
    gmds::math::Point p(1, 2, 3);
    p.setXYZ(4, 4, 4);

    REQUIRE(p.X() == Approx(4).epsilon(1e-6));
    REQUIRE(p.Y() == Approx(4).epsilon(1e-6));
    REQUIRE(p.Z() == Approx(4).epsilon(1e-6));
}

TEST_CASE("PointClass - Add") {
    gmds::math::Point p1(1, 2, 3);
    gmds::math::Point p2(2, 1, 4);

    REQUIRE((p1 + p2).X() == Approx(3).epsilon(1e-6));
    REQUIRE((p1 + p2).Y() == Approx(3).epsilon(1e-6));
    REQUIRE((p1 + p2).Z() == Approx(7).epsilon(1e-6));
}

TEST_CASE("PointClass - Minus") {
    gmds::math::Point p1(1, 2, 3);
    gmds::math::Point p2(2, 1, 3);

    REQUIRE(p1.X() - p2.X() == Approx(-1).epsilon(1e-6));
    REQUIRE(p1.Y() - p2.Y() == Approx(1).epsilon(1e-6));
    REQUIRE(p1.Z() - p2.Z()  == Approx(0).epsilon(1e-6));
}

TEST_CASE("PointClass - Multi") {
    gmds::math::Point p(1, 2, 3);
    gmds::math::Point p2 = 2 * p;
    gmds::math::Point p3 = p * 3;

    REQUIRE(p2.X() == Approx(2).epsilon(1e-6));
    REQUIRE(p2.Y() == Approx(4).epsilon(1e-6));
    REQUIRE(p2.Z() == Approx(6).epsilon(1e-6));

    REQUIRE(p3.X() == Approx(3).epsilon(1e-6));
    REQUIRE(p3.Y() == Approx(6).epsilon(1e-6));
    REQUIRE(p3.Z() == Approx(9).epsilon(1e-6));
}

TEST_CASE("PointClass - EqualAndDiff") {
    gmds::math::Point p1(2, 1, 3);
    gmds::math::Point p2(2, 1, 3);
    gmds::math::Point p3(2, 1, 4);

    REQUIRE(p1 == p2);
    REQUIRE(p1 == p1);
    REQUIRE(!(p1 != p2));
    REQUIRE(!(p1 != p1));
    REQUIRE(!(p1 < p1));
    REQUIRE(p1 < p3);
    REQUIRE(p1 <= p1);

    REQUIRE(p1.distance(p1) == Approx(0).epsilon(1e-6));
    REQUIRE(p1.distance(p2) == Approx(0).epsilon(1e-6));
    REQUIRE(p1.distance(p3) == Approx(1).epsilon(1e-6));

    REQUIRE(p1.distance2(p1) == Approx(0).epsilon(1e-6));
    REQUIRE(p1.distance2(p2) == Approx(0).epsilon(1e-6));
    REQUIRE(p1.distance2(p3) == Approx(1).epsilon(1e-6));
}

TEST_CASE("PointClass - STREAM_OUT") {
    gmds::math::Point p(2, 1, 0);
    std::stringstream s1;
    s1 << p;
    std::string s2 = "(2, 1, 0)";

    REQUIRE(s1.str() == s2);
}

TEST_CASE("PointClass - Comparison") {
    gmds::math::Point p1(2, 1, 0);
    gmds::math::Point p2(2, 2, 0);

    REQUIRE(p1 < p2);
}

TEST_CASE("PointClass - Addition") {
    gmds::math::Point p1(2, 1, -5);
    gmds::math::Point p2(2, 10, 6);

    gmds::math::Point p3 = p1 + p2;
    REQUIRE(p3.X() == p1.X() + p2.X());
    REQUIRE(p3.Y() == p1.Y() + p2.Y());
    REQUIRE(p3.Z() == p1.Z() + p2.Z());
}
