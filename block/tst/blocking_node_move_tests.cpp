#include <array>
#include <filesystem>

#include <gecko/block/Blocking.h>
#include <gecko/geom/FacetedGeometry.h>
#include <gecko/io/GmshMeshWriter.h>
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

        const auto path = (std::filesystem::temp_directory_path() / "gecko_block_node_move_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }

    /** @brief Finds the node attribute of @p ABlocking sitting exactly at @p APoint. */
    template<typename TBlocking>
    auto node_at(TBlocking &ABlocking, const Point3d &APoint) {
        auto &map = ABlocking.cmap();
        for (auto it = map.template attributes<0>().begin(), itend = map.template attributes<0>().end(); it != itend;
             ++it) {
            if (it->info().point == APoint) return it;
        }
        FAIL("no node found at the requested position");
        return map.template attributes<0>().begin();
    }
} // namespace

TEST_CASE("move_node_updates_the_node_position", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});

    const auto node = node_at(blocking, Point3d(1.0, 1.0, 0.0));
    blocking.move_node(node, Point3d(2.0, 3.0, 4.0));

    REQUIRE(node->info().point.x() == Approx(2.0));
    REQUIRE(node->info().point.y() == Approx(3.0));
    REQUIRE(node->info().point.z() == Approx(4.0));
}

TEST_CASE("move_node_refits_every_incident_edge_curve", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});

    const auto node = node_at(blocking, Point3d(1.0, 1.0, 0.0));
    blocking.move_node(node, Point3d(5.0, 1.0, 0.0));

    // Every edge touching the moved corner must now end at its new position: an edge curve left
    // pointing at the old corner would tear the block open at that corner.
    auto &map = blocking.cmap();
    int touching = 0;
    for (auto it = map.attributes<1>().begin(), itend = map.attributes<1>().end(); it != itend; ++it) {
        const auto d = it->dart();
        const auto a = map.attribute<0>(d)->info().point;
        const auto b = map.attribute<0>(map.beta<1>(d))->info().point;
        const bool touches_moved = (a == Point3d(5.0, 1.0, 0.0)) || (b == Point3d(5.0, 1.0, 0.0));
        if (!touches_moved) continue;
        ++touching;
        const auto &curve = it->info().curve;
        const Point3d start = curve.value(0.0);
        const Point3d end = curve.value(1.0);
        REQUIRE(((start == a && end == b) || (start == b && end == a)));
    }
    REQUIRE(touching == 2);
}

TEST_CASE("move_node_deforms_both_quads_sharing_the_moved_corner", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
    blocking.create_quad_block(
        {Point3d(1.0, 0.0, 0.0), Point3d(2.0, 0.0, 0.0), Point3d(2.0, 1.0, 0.0), Point3d(1.0, 1.0, 0.0)});
    blocking.build_connectivity();
    REQUIRE(blocking.nb_cells<0>() == 6); // the 2 quads really do share their common edge's corners

    const auto shared = node_at(blocking, Point3d(1.0, 1.0, 0.0));
    blocking.move_node(shared, Point3d(1.0, 4.0, 0.0));

    // Both faces' surfaces must have followed the shared corner — each evaluates to it at its own
    // local corner parameter, wherever that happens to be in its own (u,v) frame.
    auto &map = blocking.cmap();
    int faces_reaching_moved_corner = 0;
    for (auto it = map.attributes<2>().begin(), itend = map.attributes<2>().end(); it != itend; ++it) {
        const auto &surface = it->info().surface;
        for (const double u : {0.0, 1.0}) {
            for (const double v : {0.0, 1.0}) {
                if (surface.value(u, v) == Point3d(1.0, 4.0, 0.0)) {
                    ++faces_reaching_moved_corner;
                }
            }
        }
    }
    REQUIRE(faces_reaching_moved_corner == 2);
}

TEST_CASE("move_node_rebuilds_the_incident_hex_block_volume", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_hex_block({Point3d(0.0, 0.0, 0.0),
                               Point3d(1.0, 0.0, 0.0),
                               Point3d(1.0, 1.0, 0.0),
                               Point3d(0.0, 1.0, 0.0),
                               Point3d(0.0, 0.0, 1.0),
                               Point3d(1.0, 0.0, 1.0),
                               Point3d(1.0, 1.0, 1.0),
                               Point3d(0.0, 1.0, 1.0)});

    const auto node = node_at(blocking, Point3d(1.0, 1.0, 1.0));
    blocking.move_node(node, Point3d(3.0, 3.0, 3.0));

    auto &map = blocking.cmap();
    const auto block = map.attributes<3>().begin();
    bool volume_reaches_moved_corner = false;
    for (const double u : {0.0, 1.0}) {
        for (const double v : {0.0, 1.0}) {
            for (const double w : {0.0, 1.0}) {
                if (block->info().volume.value(u, v, w) == Point3d(3.0, 3.0, 3.0)) {
                    volume_reaches_moved_corner = true;
                }
            }
        }
    }
    REQUIRE(volume_reaches_moved_corner);
}

TEST_CASE("move_node_refits_curved_edges_of_a_degree_3_blocking", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});

    const auto node = node_at(blocking, Point3d(1.0, 0.0, 0.0));
    blocking.move_node(node, Point3d(4.0, 0.0, 0.0));

    // A degree-3 edge has 2 interior control points; after the move they must lie evenly along the
    // new chord, not be left where the pre-move curve had put them.
    auto &map = blocking.cmap();
    for (auto it = map.attributes<1>().begin(), itend = map.attributes<1>().end(); it != itend; ++it) {
        const auto d = it->dart();
        const auto a = map.attribute<0>(d)->info().point;
        const auto b = map.attribute<0>(map.beta<1>(d))->info().point;
        if (!(a == Point3d(0.0, 0.0, 0.0) && b == Point3d(4.0, 0.0, 0.0)) &&
            !(a == Point3d(4.0, 0.0, 0.0) && b == Point3d(0.0, 0.0, 0.0))) {
            continue;
        }
        const auto &curve = it->info().curve;
        REQUIRE(curve.value(0.5).x() == Approx(2.0));
    }
}
