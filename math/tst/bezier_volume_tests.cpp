#include <gecko/math/BezierVolume.h>
#include <gecko/math/Point3d.h>
#include <gecko/math/Vector3d.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
using Catch::Approx;
using namespace gecko;

TEST_CASE("bezier_volume_degree_1_is_trilinear", "[MathTestSuite]") {
    // Unit cube: corners [i][j][k], i/j/k along u/v/w.
    BezierVolume<1, Point3d>::ControlGrid grid;
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 2; ++j) {
            for (std::size_t k = 0; k < 2; ++k) {
                grid[i][j][k] = Point3d(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
            }
        }
    }
    const BezierVolume<1, Point3d> volume(grid);

    REQUIRE(volume.value(0.0, 0.0, 0.0) == Point3d(0.0, 0.0, 0.0));
    REQUIRE(volume.value(1.0, 1.0, 1.0) == Point3d(1.0, 1.0, 1.0));
    REQUIRE(volume.value(1.0, 0.0, 0.0) == Point3d(1.0, 0.0, 0.0));

    const Point3d center = volume.value(0.5, 0.5, 0.5);
    REQUIRE(center.x() == Approx(0.5));
    REQUIRE(center.y() == Approx(0.5));
    REQUIRE(center.z() == Approx(0.5));
}

TEST_CASE("bezier_volume_default_constructed_control_points_are_origin", "[MathTestSuite]") {
    const BezierVolume<1, Point3d> volume;
    REQUIRE(volume.control_point(0, 0, 0) == Point3d(0.0, 0.0, 0.0));
    REQUIRE(volume.control_point(1, 1, 1) == Point3d(0.0, 0.0, 0.0));
}
