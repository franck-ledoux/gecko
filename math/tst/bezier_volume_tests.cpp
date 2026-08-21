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

TEST_CASE("bezier_volume_split_along_each_axis_reproduces_the_parent_exactly", "[MathTestSuite]") {
    // The 3 axes share one implementation differing only in how a fiber maps back into the grid, so
    // this checks all 3 against an asymmetric grid — a mixed-up axis cannot survive it.
    BezierVolume<2, Point3d>::ControlGrid grid{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            for (std::size_t k = 0; k < 3; ++k) {
                grid[i][j][k] = Point3d(static_cast<double>(i) + 0.1 * static_cast<double>(k),
                                        1.5 * static_cast<double>(j) - 0.2 * static_cast<double>(i),
                                        0.7 * static_cast<double>(k) + 0.3 * static_cast<double>(i * j));
            }
        }
    }
    const BezierVolume<2, Point3d> volume(grid);

    const double s = 0.38;
    const auto [u_low, u_high] = volume.split_u(s);
    const auto [v_low, v_high] = volume.split_v(s);
    const auto [w_low, w_high] = volume.split_w(s);

    for (int a = 0; a <= 4; ++a) {
        for (int b = 0; b <= 4; ++b) {
            for (int c = 0; c <= 4; ++c) {
                const double u = static_cast<double>(a) / 4.0;
                const double v = static_cast<double>(b) / 4.0;
                const double w = static_cast<double>(c) / 4.0;
                const double hi = s + u * (1.0 - s);

                const Point3d pu = u_low.value(u, v, w);
                const Point3d eu = volume.value(u * s, v, w);
                REQUIRE(pu.x() == Approx(eu.x()).margin(1e-12));
                REQUIRE(pu.y() == Approx(eu.y()).margin(1e-12));
                REQUIRE(pu.z() == Approx(eu.z()).margin(1e-12));
                const Point3d puh = u_high.value(u, v, w);
                const Point3d euh = volume.value(hi, v, w);
                REQUIRE(puh.x() == Approx(euh.x()).margin(1e-12));
                REQUIRE(puh.z() == Approx(euh.z()).margin(1e-12));

                const Point3d pv = v_low.value(u, v, w);
                const Point3d ev = volume.value(u, v * s, w);
                REQUIRE(pv.y() == Approx(ev.y()).margin(1e-12));
                REQUIRE(pv.z() == Approx(ev.z()).margin(1e-12));

                const Point3d pw = w_low.value(u, v, w);
                const Point3d ew = volume.value(u, v, w * s);
                REQUIRE(pw.x() == Approx(ew.x()).margin(1e-12));
                REQUIRE(pw.z() == Approx(ew.z()).margin(1e-12));
            }
        }
    }
}
