#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <gmds/math/Numerics.h>
#include <gmds/math/Vector.h>

using namespace gmds;


TEST_CASE("MathTest - VectorAngle") {
    math::Vector3d v1({1, 0, 0});
    math::Vector3d v2({0, -1, 0});

    double angle1 = v1.angle(v2);
    double angle2 = v2.angle(v1);
    REQUIRE_THAT(angle1, Catch::Matchers::WithinAbs(angle2,1e-12));
}

TEST_CASE("MathTest - VectorOrientedAngle") {
    math::Vector3d v1({1, 0, 0});
    math::Vector3d v2({0, -1, 0});

    double angle1 = v1.orientedAngle(v2);
    double angle2 = v2.orientedAngle(v1);

    REQUIRE(angle1 < 0);
    REQUIRE(angle2 > 0);
    REQUIRE_THAT(angle1, Catch::Matchers::WithinAbs(-angle2,1e-12));
}

TEST_CASE("MathTest - VectorSelfAngle") {
    math::Vector3d v({1, 1, 0});

    double angle1 = v.angle(v);
    double angle2 = v.orientedAngle(v);

    REQUIRE_THAT(angle1, Catch::Matchers::WithinAbs(0.0,1e-6));
    REQUIRE_THAT(angle2, Catch::Matchers::WithinAbs(0.0,1e-6));
}

TEST_CASE("MathTest - Vector3dOrtho1") {
    math::Vector3d v({1, 0, 0});
    math::Vector3d v2 = v.getOneOrtho();

    REQUIRE(math::isZero(v.dot(v2)));
}

TEST_CASE("MathTest - Vector3dOrtho2") {
    math::Vector3d v({0, -2, 0});
    math::Vector3d v2 = v.getOneOrtho();

    REQUIRE(math::isZero(v.dot(v2)));
}

TEST_CASE("MathTest - Vector3dOrtho3") {
    math::Vector3d v({0, 0, 4});
    math::Vector3d v2 = v.getOneOrtho();

    REQUIRE(math::isZero(v.dot(v2)));
}

TEST_CASE("MathTest - Vector3dOrtho4") {
    math::Vector3d v({1.23, -0.45, 124.5});
    math::Vector3d v2 = v.getOneOrtho();

    REQUIRE(math::isZero(v.dot(v2)));
}

TEST_CASE("MathTest - VectorOrderingAngles") {
    math::Vector3d v({1, 0, 0});
    math::Vector3d v1({1, 1, 0});
    math::Vector3d v2({1, -1, 0});

    double angle1 = v.angle(v1);
    double angle2 = v.angle(v2);
    double angle1_o = v.orientedAngle(v1);
    double angle2_o = v.orientedAngle(v2);

    REQUIRE_THAT(angle2, Catch::Matchers::WithinAbs(angle1,1e-12));

    REQUIRE(angle1_o > angle2_o);
}

TEST_CASE("MathTest - numericComparison") {
    REQUIRE(math::isZero(1e-9));
    REQUIRE(math::areEquals(1, 1 + 1e-9));
}
