#include <array>
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
