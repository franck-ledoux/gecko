#include <array>
#include <filesystem>

#include <gecko/block/Blocking.h>
#include <gecko/geom/FacetedGeometry.h>
#include <gecko/io/GmshMeshWriter.h>
#include <gecko/math/Vector3d.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
using Catch::Approx;

using namespace gecko;

namespace {
    /** @brief Builds a minimal FacetedGeometry fixture: a single tagged triangle. */
    FacetedGeometry make_minimal_geom_model() {
        SimplicialMesh mesh;
        auto n0 = mesh.add_node(0, 0, 0);
        auto n1 = mesh.add_node(1, 0, 0);
        auto n2 = mesh.add_node(0, 1, 0);

        GroupRegistry groups;
        auto surf = groups.add_group("Surf", GroupDim::Dim2);
        auto &face_group = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &face_entity = mesh.add_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
        auto f0 = mesh.add_face(n0, n1, n2);
        face_group[f0.value] = surf;
        face_entity[f0.value] = 1;

        const auto path = (std::filesystem::temp_directory_path() / "gecko_block_curve_eval_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }
} // namespace

TEST_CASE("cubic_hex_block_straight_edges_evaluate_identically_to_linear", "[BlockTestSuite]") {
    // Blocking<..., BezierCurve<3,Point3d>> — degree N=3 ("curved" instantiation) — but every edge
    // is still built straight (evenly-spaced colinear control points), so it must evaluate exactly
    // like the N=1 case at every parameter, not just at the corners: this is what "N=1 collapses
    // exactly to linear blocking through the same construction" (see CoonsPatch.h) actually means
    // in practice for a real Blocking instantiation, not just the math functions in isolation.
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry, BezierCurve<3, Point3d>> blocking(geom);

    const std::array<Point3d, 8> corners = {
        Point3d(0.0, 0.0, 0.0),
        Point3d(1.0, 0.0, 0.0),
        Point3d(1.0, 1.0, 0.0),
        Point3d(0.0, 1.0, 0.0), // bottom
        Point3d(0.0, 0.0, 1.0),
        Point3d(1.0, 0.0, 1.0),
        Point3d(1.0, 1.0, 1.0),
        Point3d(0.0, 1.0, 1.0) // top
    };
    auto block = blocking.create_hex_block(corners);

    REQUIRE(blocking.nb_cells<0>() == 8);
    REQUIRE(blocking.nb_cells<1>() == 12);
    REQUIRE(blocking.nb_cells<2>() == 6);
    REQUIRE(blocking.nb_cells<3>() == 1);
    REQUIRE(blocking.is_valid_topology());

    // Same corner/center checks as the N=1 hex test — must hold exactly for N=3 too.
    REQUIRE(block->info().volume.value(0.0, 0.0, 0.0) == corners[0]);
    REQUIRE(block->info().volume.value(1.0, 1.0, 1.0) == corners[6]);
    const Point3d center = block->info().volume.value(0.5, 0.5, 0.5);
    REQUIRE(center.x() == Approx(0.5));
    REQUIRE(center.y() == Approx(0.5));
    REQUIRE(center.z() == Approx(0.5));

    // Beyond corners: a genuinely non-corner, non-center point too, to catch a curve that only
    // happens to be right at a few special parameters.
    const Point3d p = block->info().volume.value(0.25, 0.75, 0.1);
    REQUIRE(p.x() == Approx(0.25));
    REQUIRE(p.y() == Approx(0.75));
    REQUIRE(p.z() == Approx(0.1));

    // Every one of the 12 edges, individually, must trace a straight line: at every one of its
    // NumControlPoints control points (not just endpoints), the point must lie exactly on the
    // straight chord between its own first and last control point.
    for (auto it = blocking.cmap().attributes<1>().begin(), itend = blocking.cmap().attributes<1>().end(); it != itend;
         ++it) {
        const auto &curve = it->info().curve;
        const Point3d a = curve.control_points()[0];
        const Point3d b = curve.control_points()[curve.NumControlPoints - 1];
        for (std::size_t i = 0; i < curve.NumControlPoints; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(curve.NumControlPoints - 1);
            const Point3d expected = a + Vector3d(a, b) * t;
            const Point3d actual = curve.control_points()[i];
            REQUIRE(actual.x() == Approx(expected.x()).margin(1e-9));
            REQUIRE(actual.y() == Approx(expected.y()).margin(1e-9));
            REQUIRE(actual.z() == Approx(expected.z()).margin(1e-9));
        }
    }
}

TEST_CASE("cubic_quad_block_straight_surface_evaluates_identically_to_linear", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry, BezierCurve<3, Point3d>> blocking(geom);

    const std::array<Point3d, 4> corners = {
        Point3d(0.0, 0.0, 0.0), Point3d(2.0, 0.0, 0.0), Point3d(2.0, 2.0, 0.0), Point3d(0.0, 2.0, 0.0)};
    auto face = blocking.create_quad_block(corners);

    REQUIRE(blocking.nb_cells<2>() == 1);
    REQUIRE(blocking.nb_cells<3>() == 0);

    REQUIRE(face->info().surface.value(0.0, 0.0) == corners[0]);
    REQUIRE(face->info().surface.value(1.0, 1.0) == corners[2]);
    const Point3d p = face->info().surface.value(0.3, 0.7);
    REQUIRE(p.x() == Approx(0.6));
    REQUIRE(p.y() == Approx(1.4));
    REQUIRE(p.z() == Approx(0.0));
}
