#include <gecko/math/Point3d.h>
#include <gecko/math/Vector3d.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
using Catch::Approx;
using namespace gecko;

constexpr double EPSILON = 1e-12;

// =============================================================================
// Point3d TESTS
// =============================================================================

TEST_CASE("Point3d construction and element access", "[math][Point3d]") {
    SECTION("Default constructor initializes coordinates to zero") {
        constexpr Point3d p;

        REQUIRE(p.x() == Catch::Approx(0.0).margin(EPSILON));
        REQUIRE(p.y() == Catch::Approx(0.0).margin(EPSILON));
        REQUIRE(p.z() == Catch::Approx(0.0).margin(EPSILON));
    }

    SECTION("Parametrized constructor sets coordinates correctly") {
        constexpr Point3d p(1.0, 2.0, 3.0);
        REQUIRE(p.x() == Catch::Approx(1.0).margin(EPSILON));
        REQUIRE(p.y() == Catch::Approx(2.0).margin(EPSILON));
        REQUIRE(p.z() == Catch::Approx(3.0).margin(EPSILON));
    }

    SECTION("Array subscript operator provides access by index") {
        Point3d p(4.0, 5.0, 6.0);
        REQUIRE(p[0] == Catch::Approx(4.0).margin(EPSILON));
        REQUIRE(p[1] == Catch::Approx(5.0).margin(EPSILON));
        REQUIRE(p[2] == Catch::Approx(6.0).margin(EPSILON));

        p[0] = 10.0;
        p[1] = 20.0;
        p[2] = 30.0;
        REQUIRE(p.x() == Catch::Approx(10.0).margin(EPSILON));
        REQUIRE(p.y() == Catch::Approx(20.0).margin(EPSILON));
        REQUIRE(p.z() == Catch::Approx(30.0).margin(EPSILON));
    }

    SECTION("Equality comparison operator works as expected") {
        constexpr Point3d p1(1.0, 2.0, 3.0);
        constexpr Point3d p2(1.0, 2.0, 3.0);
        constexpr Point3d p3(1.0, 2.0, 3.1);

        REQUIRE(p1 == p2);
        REQUIRE_FALSE(p1 == p3);
    }
}

TEST_CASE("Affine space operations between Point3d and Vector3d", "[math][Point3d][Vector3d]") {
    SECTION("Subtracting two Point3ds yields a Vector3d") {
        const Point3d p1(4.0, 5.0, 6.0);
        const Point3d p2(1.0, 2.0, 3.0);

        const Vector3d v = p1 - p2;
        REQUIRE(v == Vector3d(3.0, 3.0, 3.0));
    }

    SECTION("Translating a Point3d by a Vector3d") {
        const Point3d p(1.0, 2.0, 3.0);
        const Vector3d v(10.0, 20.0, 30.0);

        const Point3d p_translated = p + v;
        REQUIRE(p_translated == Point3d(11.0, 22.0, 33.0));

        const Point3d p_subtracted = p - v;
        REQUIRE(p_subtracted == Point3d(-9.0, -18.0, -27.0));
    }
}

// =============================================================================
// Vector3d TESTS
// =============================================================================

TEST_CASE("Vector3d arithmetic operations", "[math][Vector3d]") {
    SECTION("Vector3d addition and subtraction") {
        const Vector3d v1(1.0, 2.0, 3.0);
        const Vector3d v2(4.0, 5.0, 6.0);

        REQUIRE((v1 + v2) == Vector3d(5.0, 7.0, 9.0));
        REQUIRE((v2 - v1) == Vector3d(3.0, 3.0, 3.0));
        REQUIRE(-v1 == Vector3d(-1.0, -2.0, -3.0));
    }

    SECTION("Scalar multiplication and division") {
        const Vector3d v(2.0, 4.0, 6.0);

        REQUIRE((v * 2.0) == Vector3d(4.0, 8.0, 12.0));
        REQUIRE((2.0 * v) == Vector3d(4.0, 8.0, 12.0));
        REQUIRE((v / 2.0) == Vector3d(1.0, 2.0, 3.0));
    }
}

TEST_CASE("Vector3d dot and cross products", "[math][Vector3d]") {
    SECTION("Dot product computation") {
        const Vector3d u(1.0, 3.0, -5.0);
        const Vector3d v(4.0, -2.0, -1.0);

        // 1*4 + 3*(-2) + (-5)*(-1) = 4 - 6 + 5 = 3
        REQUIRE(u.dot(v) == Catch::Approx(3.0).margin(EPSILON));
        REQUIRE(dot(u, v) == Catch::Approx(3.0).margin(EPSILON));
    }

    SECTION("Cross product computation and anti-commutativity") {
        const Vector3d i(1.0, 0.0, 0.0);
        const Vector3d j(0.0, 1.0, 0.0);

        // i x j = k (0, 0, 1)
        const Vector3d k = cross(i, j);
        REQUIRE(k == Vector3d(0.0, 0.0, 1.0));

        // Anti-commutativity: j x i = -k
        REQUIRE(cross(j, i) == Vector3d(0.0, 0.0, -1.0));
    }
}

TEST_CASE("Vector3d norm and normalization", "[math][Vector3d]") {
    SECTION("Calculates squared norm and norm correctly") {
        const Vector3d v(3.0, 4.0, 0.0);

        REQUIRE(v.norm_sq() == Catch::Approx(25.0).margin(EPSILON));
        REQUIRE(v.norm() == Catch::Approx(5.0).margin(EPSILON));
    }

    SECTION("Normalization yields a unit Vector3d") {
        const Vector3d v(3.0, 4.0, 0.0);
        const Vector3d unit_v = v.normalized();

        REQUIRE(unit_v.x() == Catch::Approx(0.6).margin(EPSILON));
        REQUIRE(unit_v.y() == Catch::Approx(0.8).margin(EPSILON));
        REQUIRE(unit_v.z() == Catch::Approx(0.0).margin(EPSILON));
        REQUIRE(unit_v.norm() == Catch::Approx(1.0).margin(EPSILON));
    }

    SECTION("Normalizing a zero Vector3d returns zero Vector3d safely") {
        const Vector3d zero_v(0.0, 0.0, 0.0);
        REQUIRE(zero_v.normalized() == Vector3d(0.0, 0.0, 0.0));
    }
}

TEST_CASE("Compile-time constexpr validation", "[math][constexpr]") {
    constexpr Point3d p1(1.0, 1.0, 1.0);
    constexpr Vector3d v(2.0, 0.0, 0.0);
    constexpr Point3d p2 = p1 + v;

    STATIC_REQUIRE(p2.x() == 3.0);
    STATIC_REQUIRE(dot(v, v) == 4.0);
}