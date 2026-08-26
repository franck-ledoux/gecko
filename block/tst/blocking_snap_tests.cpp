#include <array>
#include <cmath>
#include <numbers>
#include <filesystem>
#include <utility>
#include <vector>

#include <gecko/block/Blocking.h>
#include <gecko/geom/FacetedGeometry.h>
#include <gecko/io/GmshMeshWriter.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
using Catch::Approx;

using namespace gecko;

namespace {
    /**
     * @brief The same unit-square boundary rep as blocking_classification_tests.cpp's own fixture:
     * 4 tagged vertices (1-4) at the corners, 4 tagged straight boundary curves (10 bottom, 11
     * right, 12 top, 13 left) and 1 tagged surface (20), all in the z=0 plane.
     *
     * Duplicated rather than shared because these tests exercise a different property — that a
     * cell's classification follows its *boundary's* classification — and want to state the
     * incidence they rely on (which vertex bounds which curve) right where it is asserted.
     */
    FacetedGeometry make_square_geom_model() {
        SimplicialMesh mesh;
        auto v0 = mesh.add_node(0, 0, 0);
        auto v1 = mesh.add_node(1, 0, 0);
        auto v2 = mesh.add_node(1, 1, 0);
        auto v3 = mesh.add_node(0, 1, 0);

        GroupRegistry groups;
        auto vtx_group = groups.add_group("Vertices", GroupDim::Dim0);
        auto curve_group = groups.add_group("Curves", GroupDim::Dim1);
        auto surf_group = groups.add_group("Surf", GroupDim::Dim2);

        auto &node_group = mesh.add_variable<GroupId, CellType::Node>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &node_entity = mesh.add_variable<Int, CellType::Node>(std::string(io::ENTITY_TAG_VARIABLE));
        const std::array<NodeId, 4> corners{v0, v1, v2, v3};
        for (std::size_t i = 0; i < 4; ++i) {
            node_group[corners[i].value] = vtx_group;
            node_entity[corners[i].value] = static_cast<Int>(i + 1);
        }

        auto &edge_group = mesh.add_variable<GroupId, CellType::Edge>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &edge_entity = mesh.add_variable<Int, CellType::Edge>(std::string(io::ENTITY_TAG_VARIABLE));
        const std::array<std::pair<NodeId, NodeId>, 4> curve_edges{
            std::pair{v0, v1}, std::pair{v1, v2}, std::pair{v2, v3}, std::pair{v3, v0}};
        for (std::size_t i = 0; i < 4; ++i) {
            auto e = mesh.add_edge(curve_edges[i].first, curve_edges[i].second);
            edge_group[e.value] = curve_group;
            edge_entity[e.value] = static_cast<Int>(10 + i);
        }

        auto &face_group = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &face_entity = mesh.add_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
        for (auto f : {mesh.add_face(v0, v1, v2), mesh.add_face(v0, v2, v3)}) {
            face_group[f.value] = surf_group;
            face_entity[f.value] = 20;
        }

        const auto path = (std::filesystem::temp_directory_path() / "gecko_block_snap_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }

    /** @brief The single edge of @p ABlocking joining the corners at @p AP0 and @p AP1. */
    template<typename TBlocking>
    typename TBlocking::Edge edge_between(TBlocking &ABlocking, const Point3d &AP0, const Point3d &AP1) {
        for (auto it = ABlocking.cmap().template attributes<1>().begin(),
                  itend = ABlocking.cmap().template attributes<1>().end();
             it != itend;
             ++it) {
            const auto d = it->dart();
            const Point3d &a = ABlocking.cmap().template attribute<0>(d)->info().point;
            const Point3d &b =
                ABlocking.cmap().template attribute<0>(ABlocking.cmap().template beta<1>(d))->info().point;
            if ((a == AP0 && b == AP1) || (a == AP1 && b == AP0)) return it;
        }
        FAIL("no edge joins those two corners");
        return {};
    }

    /** @brief The node of @p ABlocking sitting at @p AP. */
    template<typename TBlocking>
    typename TBlocking::Node node_at(TBlocking &ABlocking, const Point3d &AP) {
        for (auto it = ABlocking.cmap().template attributes<0>().begin(),
                  itend = ABlocking.cmap().template attributes<0>().end();
             it != itend;
             ++it) {
            if (it->info().point == AP) return it;
        }
        FAIL("no node at that position");
        return {};
    }
} // namespace

TEST_CASE("edge_between_two_vertices_of_one_curve_is_classified_on_that_curve", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_quad_block({Point3d(0, 0, 0), Point3d(1, 0, 0), Point3d(1, 1, 0), Point3d(0, 1, 0)});
    blocking.classify(1e-6);

    // Corners (0,0,0) and (1,0,0) are vertices 1 and 2, which share exactly one curve: the bottom
    // one, tag 10. That is the lowest-dimensional entity containing both.
    const auto e = edge_between(blocking, Point3d(0, 0, 0), Point3d(1, 0, 0));
    REQUIRE(e->info().geom_targets.size() == 1);
    REQUIRE(e->info().geom_targets[0] == std::pair{GroupDim::Dim1, Int{10}});
}

TEST_CASE("diagonal_edge_is_classified_on_the_surface_not_a_boundary_curve", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    // Corners 0 and 1 are consecutive, so the block edge between them runs (0,0,0)->(1,1,0): right
    // across the square's interior, along no boundary curve. The other 2 corners are far enough
    // away to stay unclassified, which leaves this edge the only one inference can decide.
    blocking.create_quad_block({Point3d(0, 0, 0), Point3d(1, 1, 0), Point3d(2, 2, 0), Point3d(-1, 1, 0)});

    // A tolerance loose enough to reach the surface also puts all 4 boundary curves within reach of
    // the diagonal's midpoint (0.5,0.5,0) — each is exactly 0.5 away. Classifying by that midpoint
    // alone would therefore stop at dimension 1 and pin the diagonal to a boundary curve it does
    // not lie on. Its 2 corners are vertices 1 and 3, whose only common entity is surface 20.
    blocking.classify(1e-6, 0.6, 0.6);

    const auto e = edge_between(blocking, Point3d(0, 0, 0), Point3d(1, 1, 0));
    REQUIRE(e->info().geom_targets.size() == 1);
    REQUIRE(e->info().geom_targets[0] == std::pair{GroupDim::Dim2, Int{20}});
}

TEST_CASE("face_is_classified_on_the_surface_containing_all_its_edges", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_quad_block({Point3d(0, 0, 0), Point3d(1, 0, 0), Point3d(1, 1, 0), Point3d(0, 1, 0)});
    blocking.classify(1e-6);

    // The 4 edges land on the 4 boundary curves; the only entity containing all of them is the
    // surface.
    auto f = blocking.cmap().template attributes<2>().begin();
    REQUIRE(f->info().geom_targets.size() == 1);
    REQUIRE(f->info().geom_targets[0] == std::pair{GroupDim::Dim2, Int{20}});
}

TEST_CASE("edge_with_an_unclassified_corner_falls_back_to_proximity", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    // The far corner is nowhere near the model, so it stays unclassified and can constrain nothing.
    blocking.create_quad_block({Point3d(0, 0, 0), Point3d(1, 0, 0), Point3d(9, 9, 9), Point3d(0, 1, 0)});
    blocking.classify(1e-6);

    const auto e = edge_between(blocking, Point3d(1, 0, 0), Point3d(9, 9, 9));
    REQUIRE(e->info().geom_targets.empty()); // Nothing within 1e-6 of its midpoint either.

    // The edge between the 2 classified corners is unaffected by its neighbour's failure.
    const auto ok = edge_between(blocking, Point3d(0, 0, 0), Point3d(1, 0, 0));
    REQUIRE(ok->info().geom_targets.size() == 1);
    REQUIRE(ok->info().geom_targets[0] == std::pair{GroupDim::Dim1, Int{10}});
}

TEST_CASE("snap_node_reclassifies_the_node_and_its_incident_edges", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_quad_block({Point3d(0, 0, 0), Point3d(1, 0, 0), Point3d(1, 1, 0), Point3d(0, 1, 0)});
    blocking.classify(1e-6);

    // Drag the (1,1,0) corner off the model, into the middle of the bottom curve's neighbourhood.
    auto n = node_at(blocking, Point3d(1, 1, 0));
    blocking.move_node(n, Point3d(0.5, 0.02, 0));
    blocking.snap_node(n, 0.1, 0.1, 0.1);

    // It snapped onto the bottom curve (tag 10), not onto a vertex: (0.5,0.02,0) is 0.5 away from
    // the nearest corner but 0.02 from the curve.
    REQUIRE(n->info().geom_targets.size() == 1);
    REQUIRE(n->info().geom_targets[0] == std::pair{GroupDim::Dim1, Int{10}});
    REQUIRE(n->info().point.y() == Approx(0.0));

    // Its edge to (1,0,0) — vertex 2, which also bounds curve 10 — now lies along that same curve.
    const auto e = edge_between(blocking, Point3d(0.5, 0, 0), Point3d(1, 0, 0));
    REQUIRE(e->info().geom_targets.size() == 1);
    REQUIRE(e->info().geom_targets[0] == std::pair{GroupDim::Dim1, Int{10}});
}

TEST_CASE("project_onto_classification_reads_where_a_node_already_is_not_where_it_is_near", "[BlockTestSuite]") {
    // The read-only counterpart of snap_node()'s own search: it answers from a node's *current*
    // classification, and must not re-run the search itself — which the fixture checks by asking
    // from a point that is actually closer to a different vertex than the one the node is on.
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_quad_block({Point3d(0, 0, 0), Point3d(1, 0, 0), Point3d(1, 1, 0), Point3d(0, 1, 0)});
    blocking.classify(1e-6);

    const auto origin = node_at(blocking, Point3d(0, 0, 0));
    REQUIRE(origin->info().geom_targets.size() == 1);
    REQUIRE(origin->info().geom_targets[0] == std::pair{GroupDim::Dim0, Int{1}});

    // Nearer to (1,1,0) than to the origin, yet the answer is still the origin: it reads what
    // `origin` is classified on, not what the trial point happens to be near.
    const Point3d trial(0.9, 0.9, 5.0);
    const Point3d projected = blocking.project_onto_classification(origin, trial);
    REQUIRE(projected.x() == Approx(0.0));
    REQUIRE(projected.y() == Approx(0.0));
    REQUIRE(projected.z() == Approx(0.0));

    // An unclassified node gives the trial point straight back — nothing to pull it onto.
    const auto bare = blocking.create_node(Point3d(3, 3, 3));
    REQUIRE(bare->info().geom_targets.empty());
    const Point3d unmoved = blocking.project_onto_classification(bare, trial);
    REQUIRE(unmoved.x() == Approx(trial.x()));
    REQUIRE(unmoved.y() == Approx(trial.y()));
    REQUIRE(unmoved.z() == Approx(trial.z()));
}

TEST_CASE("a_curved_edge_lies_on_its_geometric_curve_not_merely_near_it", "[BlockTestSuite]") {
    // The bottom boundary curve is a 2-segment polyline bulging up to (0.5, 0.25, 0), so following
    // it genuinely requires bending — a straight edge is 0.25 away from it at mid-span.
    SimplicialMesh mesh;
    auto v0 = mesh.add_node(0, 0, 0);
    auto v1 = mesh.add_node(1, 0, 0);
    auto v2 = mesh.add_node(1, 1, 0);
    auto v3 = mesh.add_node(0, 1, 0);
    auto vm = mesh.add_node(0.5, 0.25, 0);

    GroupRegistry groups;
    auto vtx_group = groups.add_group("Vertices", GroupDim::Dim0);
    auto curve_group = groups.add_group("Curves", GroupDim::Dim1);
    auto surf_group = groups.add_group("Surf", GroupDim::Dim2);

    auto &node_group = mesh.add_variable<GroupId, CellType::Node>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    auto &node_entity = mesh.add_variable<Int, CellType::Node>(std::string(io::ENTITY_TAG_VARIABLE));
    const std::array<NodeId, 4> corners{v0, v1, v2, v3};
    for (std::size_t i = 0; i < 4; ++i) {
        node_group[corners[i].value] = vtx_group;
        node_entity[corners[i].value] = static_cast<Int>(i + 1);
    }

    auto &edge_group = mesh.add_variable<GroupId, CellType::Edge>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    auto &edge_entity = mesh.add_variable<Int, CellType::Edge>(std::string(io::ENTITY_TAG_VARIABLE));
    for (auto e : {mesh.add_edge(v0, vm), mesh.add_edge(vm, v1)}) {
        edge_group[e.value] = curve_group;
        edge_entity[e.value] = 10;
    }

    auto &face_group = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    auto &face_entity = mesh.add_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
    for (auto f : {mesh.add_face(v0, v1, v2), mesh.add_face(v0, v2, v3)}) {
        face_group[f.value] = surf_group;
        face_entity[f.value] = 20;
    }

    const auto path = (std::filesystem::temp_directory_path() / "gecko_block_bent_curve_test.msh").string();
    io::SimplicialMeshWriter::write(path, mesh, groups);
    const FacetedGeometry geom(path);
    std::filesystem::remove(path);

    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_quad_block({Point3d(0, 0, 0), Point3d(1, 0, 0), Point3d(1, 1, 0), Point3d(0, 1, 0)});
    blocking.classify(1e-6);

    const auto e = edge_between(blocking, Point3d(0, 0, 0), Point3d(1, 0, 0));
    REQUIRE(e->info().geom_targets[0] == std::pair{GroupDim::Dim1, Int{10}});

    const auto *curve = geom.curve_by_tag(10);
    REQUIRE(curve != nullptr);

    // No smooth cubic reproduces this polyline's kink, so the fit can only approximate it. Measure
    // that residual against what projecting the control points would have given on the very same
    // edge, rather than trusting a hand-picked threshold.
    const auto worst_distance = [&curve](const BezierCurve<Point3d> &ACurve) {
        double worst = 0.0;
        for (int i = 0; i <= 50; ++i) {
            worst = std::max(worst, curve->distance(ACurve.value(static_cast<double>(i) / 50.0)));
        }
        return worst;
    };

    BezierCurve<Point3d> control_point_projection(3);
    const Vector3d chord(Point3d(0, 0, 0), Point3d(1, 0, 0));
    for (std::size_t i = 0; i < 4; ++i) {
        Point3d cp = Point3d(0, 0, 0) + chord * (static_cast<double>(i) / 3.0);
        if (i > 0 && i < 3) curve->project(cp);
        control_point_projection[i] = cp;
    }

    // Measured here: 0.089 interpolated vs 0.134 projected. The gap is modest on this deliberately
    // pathological V — its kink is unreachable by any smooth cubic — and the residual is a sampling
    // limitation, not a fitting one: evenly spaced chord points never sample near the apex. On the
    // smooth curves real CAD produces, the interpolated fit is essentially exact. What matters is
    // that the curve now touches the geometry at every sample parameter, which the assertions above
    // pin exactly, instead of only at its 2 endpoints.
    const double fitted_error = worst_distance(e->info().curve);
    const double projected_error = worst_distance(control_point_projection);
    INFO("interpolated: " << fitted_error << "  control-point projection: " << projected_error);
    REQUIRE(fitted_error < projected_error);

    // And it really is bent: the straight chord it started as passes 0.25 from the apex.
    REQUIRE(e->info().curve.value(0.5).y() > 0.14);
}

TEST_CASE("a_curved_edge_leaves_its_ends_along_the_geometric_curve", "[BlockTestSuite]") {
    // A quarter circle of radius 1, finely faceted: smooth, like real CAD, and with a tangent known
    // in closed form — at (cos a, sin a) it is (-sin a, cos a).
    constexpr int SEGMENTS = 64;
    SimplicialMesh mesh;
    std::vector<NodeId> arc;
    for (int i = 0; i <= SEGMENTS; ++i) {
        const double a = (std::numbers::pi / 2.0) * static_cast<double>(i) / SEGMENTS;
        arc.push_back(mesh.add_node(std::cos(a), std::sin(a), 0));
    }
    const auto apex = mesh.add_node(2, 2, 0);

    GroupRegistry groups;
    const auto vtx_group = groups.add_group("Vertices", GroupDim::Dim0);
    const auto curve_group = groups.add_group("Curves", GroupDim::Dim1);
    const auto surf_group = groups.add_group("Surf", GroupDim::Dim2);

    auto &node_group = mesh.add_variable<GroupId, CellType::Node>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    auto &node_entity = mesh.add_variable<Int, CellType::Node>(std::string(io::ENTITY_TAG_VARIABLE));
    node_group[arc.front().value] = vtx_group;
    node_entity[arc.front().value] = 1;
    node_group[arc.back().value] = vtx_group;
    node_entity[arc.back().value] = 2;

    auto &edge_group = mesh.add_variable<GroupId, CellType::Edge>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    auto &edge_entity = mesh.add_variable<Int, CellType::Edge>(std::string(io::ENTITY_TAG_VARIABLE));
    auto &face_group = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    auto &face_entity = mesh.add_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
    for (int i = 0; i < SEGMENTS; ++i) {
        auto e = mesh.add_edge(arc[i], arc[i + 1]);
        edge_group[e.value] = curve_group;
        edge_entity[e.value] = 10;
        auto f = mesh.add_face(arc[i], arc[i + 1], apex);
        face_group[f.value] = surf_group;
        face_entity[f.value] = 20;
    }

    const auto path = (std::filesystem::temp_directory_path() / "gecko_block_arc_test.msh").string();
    io::SimplicialMeshWriter::write(path, mesh, groups);
    const FacetedGeometry geom(path);
    std::filesystem::remove(path);

    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_quad_block({Point3d(1, 0, 0), Point3d(0, 1, 0), Point3d(2, 2, 0), Point3d(2, 0, 0)});
    // Grabbed before classifying, since snapping nudges the corner onto the arc's own end node —
    // which sits at (cos(pi/2), sin(pi/2)), not exactly (0,1,0), and so no longer matches by value.
    const auto e = edge_between(blocking, Point3d(1, 0, 0), Point3d(0, 1, 0));
    blocking.classify(1e-6, 0.5, 0.5);

    REQUIRE(e->info().geom_targets[0] == std::pair{GroupDim::Dim1, Int{10}});

    // A Bezier leaves its start along P1-P0. Interpolating positions alone leaves that direction
    // unconstrained, and it comes out ~30 degrees off here; the fit pins it to the geometry's own
    // tangent instead. What is left is the faceting's angular resolution: 90 degrees over 64
    // segments means a chord direction sits ~0.7 degrees off the true tangent.
    const auto &cp = e->info().curve.control_points();
    const Vector3d start = Vector3d(cp[0], cp[1]).normalized();
    const Vector3d end = Vector3d(cp[2], cp[3]).normalized();
    REQUIRE(start.dot(Vector3d(0, 1, 0)) > std::cos(1.5 * std::numbers::pi / 180.0));
    REQUIRE(end.dot(Vector3d(-1, 0, 0)) > std::cos(1.5 * std::numbers::pi / 180.0));

    // Getting the ends right is most of getting the shape right: the worst deviation over the whole
    // span measures 0.006 on this unit-radius arc, against 0.034 for a fit free to pick its own end
    // tangents.
    const auto *curve = geom.curve_by_tag(10);
    for (int i = 0; i <= 100; ++i) {
        REQUIRE(curve->distance(e->info().curve.value(static_cast<double>(i) / 100.0)) < 0.01);
    }
}

TEST_CASE("snap_node_leaves_cells_it_does_not_touch_alone", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_quad_block({Point3d(0, 0, 0), Point3d(1, 0, 0), Point3d(1, 1, 0), Point3d(0, 1, 0)});
    blocking.classify(1e-6);

    // The top edge (1,1,0)-(0,1,0) is on curve 12 and touches neither of the corners below.
    const auto untouched = edge_between(blocking, Point3d(1, 1, 0), Point3d(0, 1, 0));
    const auto before = untouched->info().geom_targets;

    auto n = node_at(blocking, Point3d(0, 0, 0));
    blocking.move_node(n, Point3d(0.02, 0.02, 0));
    blocking.snap_node(n, 0.1, 0.1, 0.1);

    REQUIRE(untouched->info().geom_targets == before);
}
