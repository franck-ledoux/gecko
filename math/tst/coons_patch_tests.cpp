#include <gecko/math/BezierCurve.h>
#include <gecko/math/BezierSurface.h>
#include <gecko/math/BezierVolume.h>
#include <gecko/math/CoonsPatch.h>
#include <gecko/math/Point3d.h>
#include <gecko/math/Vector3d.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
using Catch::Approx;
using namespace gecko;

TEST_CASE("coons_surface_from_edges_degree_1_is_bilinear", "[MathTestSuite]") {
    const Point3d p00(0.0, 0.0, 0.0), p10(1.0, 0.0, 0.0), p01(0.0, 1.0, 0.0), p11(1.0, 1.0, 0.0);
    const BezierCurve<Point3d> edge_u0(p00, p10); // v=0
    const BezierCurve<Point3d> edge_u1(p01, p11); // v=1
    const BezierCurve<Point3d> edge_v0(p00, p01); // u=0
    const BezierCurve<Point3d> edge_v1(p10, p11); // u=1

    const auto surface = coons_surface_from_edges(edge_u0, edge_u1, edge_v0, edge_v1);

    REQUIRE(surface.control_point(0, 0) == p00);
    REQUIRE(surface.control_point(1, 0) == p10);
    REQUIRE(surface.control_point(0, 1) == p01);
    REQUIRE(surface.control_point(1, 1) == p11);

    const Point3d center = surface.value(0.5, 0.5);
    REQUIRE(center.x() == Approx(0.5));
    REQUIRE(center.y() == Approx(0.5));
    REQUIRE(center.z() == Approx(0.0));
}

TEST_CASE("coons_surface_from_edges_reproduces_a_bulging_boundary", "[MathTestSuite]") {
    const Point3d p00(0.0, 0.0, 0.0), p10(2.0, 0.0, 0.0), p01(0.0, 1.0, 0.0), p11(2.0, 1.0, 0.0);
    // v=0 boundary bulges up in z; the 3 other boundaries stay straight (degree-2 for consistency).
    const BezierCurve<Point3d> edge_u0(p00, Point3d(1.0, 0.0, 3.0), p10);
    const BezierCurve<Point3d> edge_u1(p01, Point3d(1.0, 1.0, 0.0), p11);
    const BezierCurve<Point3d> edge_v0(p00, Point3d(0.0, 0.5, 0.0), p01);
    const BezierCurve<Point3d> edge_v1(p10, Point3d(2.0, 0.5, 0.0), p11);

    const auto surface = coons_surface_from_edges(edge_u0, edge_u1, edge_v0, edge_v1);

    // The surface reproduces edge_u0 exactly along its v=0 boundary (j=0).
    for (std::size_t i = 0; i < 3; ++i) {
        REQUIRE(surface.control_point(i, 0) == edge_u0.control_points()[i]);
    }
    // Near v=0 the surface must still bulge up in z, following edge_u0's own bulge.
    REQUIRE(surface.control_point(1, 0).z() > 1.0);
}

TEST_CASE("tfi_volume_from_faces_degree_1_is_trilinear", "[MathTestSuite]") {
    // 8 corners of a unit cube, c[i][j][k].
    Point3d c[2][2][2];
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
                c[i][j][k] = Point3d(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));

    // 12 edges, one per (axis, pair of fixed indices on the other 2 axes).
    const BezierCurve<Point3d> ei00(c[0][0][0], c[1][0][0]), ei01(c[0][0][1], c[1][0][1]), ei10(c[0][1][0], c[1][1][0]),
        ei11(c[0][1][1], c[1][1][1]);
    const BezierCurve<Point3d> ej00(c[0][0][0], c[0][1][0]), ej01(c[0][0][1], c[0][1][1]), ej10(c[1][0][0], c[1][1][0]),
        ej11(c[1][0][1], c[1][1][1]);
    const BezierCurve<Point3d> ek00(c[0][0][0], c[0][0][1]), ek01(c[0][1][0], c[0][1][1]), ek10(c[1][0][0], c[1][0][1]),
        ek11(c[1][1][0], c[1][1][1]);

    // 6 bounding faces, each a Coons patch of 4 of the 12 edges (grid indexing per CoonsPatch.h's doc).
    const auto face_u0 = coons_surface_from_edges(ej00, ej01, ek00, ek01); // u=0, grid[v][w]
    const auto face_u1 = coons_surface_from_edges(ej10, ej11, ek10, ek11); // u=1, grid[v][w]
    const auto face_v0 = coons_surface_from_edges(ei00, ei01, ek00, ek10); // v=0, grid[u][w]
    const auto face_v1 = coons_surface_from_edges(ei10, ei11, ek01, ek11); // v=1, grid[u][w]
    const auto face_w0 = coons_surface_from_edges(ei00, ei10, ej00, ej10); // w=0, grid[u][v]
    const auto face_w1 = coons_surface_from_edges(ei01, ei11, ej01, ej11); // w=1, grid[u][v]

    const auto volume = tfi_volume_from_faces(face_u0, face_u1, face_v0, face_v1, face_w0, face_w1);

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                REQUIRE(volume.control_point(static_cast<std::size_t>(i),
                                             static_cast<std::size_t>(j),
                                             static_cast<std::size_t>(k)) == c[i][j][k]);
            }
        }
    }

    const Point3d center = volume.value(0.5, 0.5, 0.5);
    REQUIRE(center.x() == Approx(0.5));
    REQUIRE(center.y() == Approx(0.5));
    REQUIRE(center.z() == Approx(0.5));
}
