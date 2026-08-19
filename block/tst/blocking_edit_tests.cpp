#include <array>
#include <filesystem>

#include <gecko/block/Blocking.h>
#include <gecko/geom/FacetedGeometry.h>
#include <gecko/io/GmshMeshWriter.h>
#include <catch2/catch_test_macros.hpp>

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

        const auto path = (std::filesystem::temp_directory_path() / "gecko_block_edit_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }
} // namespace

TEST_CASE("is_purely_2d_true_for_quad_only_blocking", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
    REQUIRE(blocking.is_purely_2d());
}

TEST_CASE("is_purely_2d_false_once_a_hex_block_exists", "[BlockTestSuite]") {
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
    REQUIRE_FALSE(blocking.is_purely_2d());

    // can_delete_face must reject too, even given a legitimate face handle (one of the hex's own
    // bounding faces) — the 2D-only precondition applies regardless of which face is targeted.
    auto face_it = blocking.cmap().attributes<2>().begin();
    REQUIRE(face_it != blocking.cmap().attributes<2>().end());
    REQUIRE_FALSE(blocking.can_delete_face(face_it));
}

TEST_CASE("delete_face_removes_face_and_garbage_collects_its_unshared_boundary", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    // Quad A: [0,1]x[0,1]. Quad B: [1,2]x[0,1], sharing A's edge (1,0)-(1,1) with B's edge (1,1)-(1,0).
    const std::array<Point3d, 4> corners_a = {
        Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)};
    const std::array<Point3d, 4> corners_b = {
        Point3d(1.0, 0.0, 0.0), Point3d(2.0, 0.0, 0.0), Point3d(2.0, 1.0, 0.0), Point3d(1.0, 1.0, 0.0)};

    const auto face_a = blocking.create_quad_block(corners_a);
    blocking.create_quad_block(corners_b);
    blocking.build_connectivity();

    // Baseline, matching blocking_connectivity_tests.cpp's equivalent fixture: 6 nodes, 7 edges
    // (1 shared), 2 faces.
    REQUIRE(blocking.nb_cells<0>() == 6);
    REQUIRE(blocking.nb_cells<1>() == 7);
    REQUIRE(blocking.nb_cells<2>() == 2);

    REQUIRE(blocking.can_delete_face(face_a));
    blocking.delete_face(face_a);

    REQUIRE(blocking.is_valid_topology());
    // Face A is gone (2-1=1 face left: B). Of A's 4 boundary edges, exactly 1 (the shared one) is
    // still referenced by B and survives; the other 3 (unshared) are garbage-collected by CGAL's
    // own attribute reference counting (no dart left referencing them) -> 7-3=4 edges. Of A's 4
    // corners, exactly 2 (the shared edge's endpoints) still belong to B and survive; the other 2
    // are garbage-collected -> 6-2=4 nodes.
    REQUIRE(blocking.nb_cells<0>() == 4);
    REQUIRE(blocking.nb_cells<1>() == 4);
    REQUIRE(blocking.nb_cells<2>() == 1);

    // The survivor really is B, intact: exactly 1 face remains, and every one of B's 4 corners is
    // present among the 4 surviving nodes.
    REQUIRE(blocking.cmap().attributes<2>().begin() != blocking.cmap().attributes<2>().end());
    for (const Point3d &expected : corners_b) {
        bool present = false;
        for (auto it = blocking.cmap().attributes<0>().begin(), itend = blocking.cmap().attributes<0>().end();
             it != itend;
             ++it) {
            if (it->info().point == expected) {
                present = true;
                break;
            }
        }
        REQUIRE(present);
    }
}

TEST_CASE("delete_face_on_a_standalone_face_removes_everything", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    const auto face = blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});

    REQUIRE(blocking.nb_cells<0>() == 4);
    REQUIRE(blocking.nb_cells<1>() == 4);
    REQUIRE(blocking.nb_cells<2>() == 1);

    REQUIRE(blocking.can_delete_face(face));
    blocking.delete_face(face);

    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.nb_cells<0>() == 0);
    REQUIRE(blocking.nb_cells<1>() == 0);
    REQUIRE(blocking.nb_cells<2>() == 0);
}
