#include <algorithm>
#include <array>
#include <filesystem>
#include <vector>

#include <unit_test_config.h>

#include <gecko/block/Blocking.h>
#include <gecko/geom/FacetedGeometry.h>
#include <gecko/io/GmshMeshWriter.h>
#include <gecko/math/Vector3d.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
using Catch::Approx;

using namespace gecko;

namespace {
    /**
     * @brief A unit-square boundary-rep fixture: 4 tagged vertices (1-4), 4 tagged straight
     * boundary curves (10-13, one edge each) and 1 tagged surface (20, 2 triangles), all in the
     * z=0 plane, corners in perimeter order (0,0,0)->(1,0,0)->(1,1,0)->(0,1,0).
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
        node_group[v0.value] = vtx_group;
        node_entity[v0.value] = 1;
        node_group[v1.value] = vtx_group;
        node_entity[v1.value] = 2;
        node_group[v2.value] = vtx_group;
        node_entity[v2.value] = 3;
        node_group[v3.value] = vtx_group;
        node_entity[v3.value] = 4;

        auto &edge_group = mesh.add_variable<GroupId, CellType::Edge>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &edge_entity = mesh.add_variable<Int, CellType::Edge>(std::string(io::ENTITY_TAG_VARIABLE));
        auto e_bottom = mesh.add_edge(v0, v1);
        edge_group[e_bottom.value] = curve_group;
        edge_entity[e_bottom.value] = 10;
        auto e_right = mesh.add_edge(v1, v2);
        edge_group[e_right.value] = curve_group;
        edge_entity[e_right.value] = 11;
        auto e_top = mesh.add_edge(v2, v3);
        edge_group[e_top.value] = curve_group;
        edge_entity[e_top.value] = 12;
        auto e_left = mesh.add_edge(v3, v0);
        edge_group[e_left.value] = curve_group;
        edge_entity[e_left.value] = 13;

        auto &face_group = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &face_entity = mesh.add_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
        auto f0 = mesh.add_face(v0, v1, v2);
        face_group[f0.value] = surf_group;
        face_entity[f0.value] = 20;
        auto f1 = mesh.add_face(v0, v2, v3);
        face_group[f1.value] = surf_group;
        face_entity[f1.value] = 20;

        const auto path = (std::filesystem::temp_directory_path() / "gecko_block_classify_square_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }

    /**
     * @brief Same square as make_square_geom_model(), except the bottom boundary curve (tag 10) is
     * a 2-segment polyline bulging to (0.5, 0.4, 0) instead of a straight segment — used to verify
     * that classify() actually reshapes a block edge (and, through it, its owning face) to follow a
     * genuinely curved reference curve, not just snap onto a straight one.
     */
    FacetedGeometry make_square_geom_model_with_bent_bottom() {
        SimplicialMesh mesh;
        auto v0 = mesh.add_node(0, 0, 0);
        auto v1 = mesh.add_node(1, 0, 0);
        auto v2 = mesh.add_node(1, 1, 0);
        auto v3 = mesh.add_node(0, 1, 0);
        auto vm = mesh.add_node(0.5, 0.4, 0);

        GroupRegistry groups;
        auto vtx_group = groups.add_group("Vertices", GroupDim::Dim0);
        auto curve_group = groups.add_group("Curves", GroupDim::Dim1);
        auto surf_group = groups.add_group("Surf", GroupDim::Dim2);

        auto &node_group = mesh.add_variable<GroupId, CellType::Node>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &node_entity = mesh.add_variable<Int, CellType::Node>(std::string(io::ENTITY_TAG_VARIABLE));
        node_group[v0.value] = vtx_group;
        node_entity[v0.value] = 1;
        node_group[v1.value] = vtx_group;
        node_entity[v1.value] = 2;
        node_group[v2.value] = vtx_group;
        node_entity[v2.value] = 3;
        node_group[v3.value] = vtx_group;
        node_entity[v3.value] = 4;

        auto &edge_group = mesh.add_variable<GroupId, CellType::Edge>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &edge_entity = mesh.add_variable<Int, CellType::Edge>(std::string(io::ENTITY_TAG_VARIABLE));
        auto e_bottom0 = mesh.add_edge(v0, vm);
        edge_group[e_bottom0.value] = curve_group;
        edge_entity[e_bottom0.value] = 10;
        auto e_bottom1 = mesh.add_edge(vm, v1);
        edge_group[e_bottom1.value] = curve_group;
        edge_entity[e_bottom1.value] = 10;
        auto e_right = mesh.add_edge(v1, v2);
        edge_group[e_right.value] = curve_group;
        edge_entity[e_right.value] = 11;
        auto e_top = mesh.add_edge(v2, v3);
        edge_group[e_top.value] = curve_group;
        edge_entity[e_top.value] = 12;
        auto e_left = mesh.add_edge(v3, v0);
        edge_group[e_left.value] = curve_group;
        edge_entity[e_left.value] = 13;

        auto &face_group = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &face_entity = mesh.add_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
        auto f0 = mesh.add_face(v0, v1, v2);
        face_group[f0.value] = surf_group;
        face_entity[f0.value] = 20;
        auto f1 = mesh.add_face(v0, v2, v3);
        face_group[f1.value] = surf_group;
        face_entity[f1.value] = 20;

        const auto path =
            (std::filesystem::temp_directory_path() / "gecko_block_classify_bent_square_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }

    /**
     * @brief A geometric model with no vertices at all: 2 straight curves (100, 101) sharing an
     * untagged endpoint at the origin — used to exercise the "classified on a *group* of entities"
     * design (see CellInfo::geom_targets doc): a block node placed near that untagged junction has
     * no vertex to snap to, so its search falls through to dimension 1 and must find *both* curves.
     */
    FacetedGeometry make_curve_junction_geom_model() {
        SimplicialMesh mesh;
        auto n0 = mesh.add_node(0, 0, 0);
        auto n1 = mesh.add_node(1, 0, 0);
        auto n2 = mesh.add_node(0, 1, 0);

        GroupRegistry groups;
        auto curve_group = groups.add_group("Curves", GroupDim::Dim1);

        auto &edge_group = mesh.add_variable<GroupId, CellType::Edge>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &edge_entity = mesh.add_variable<Int, CellType::Edge>(std::string(io::ENTITY_TAG_VARIABLE));
        auto e0 = mesh.add_edge(n0, n1);
        edge_group[e0.value] = curve_group;
        edge_entity[e0.value] = 100;
        auto e1 = mesh.add_edge(n0, n2);
        edge_group[e1.value] = curve_group;
        edge_entity[e1.value] = 101;

        const auto path = (std::filesystem::temp_directory_path() / "gecko_block_classify_junction_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }

    /** @brief Checks whether @p AP lies within @p ATol of any point in @p ACandidates. */
    bool is_near_any(const Point3d &AP, const std::vector<Point3d> &ACandidates, double ATol) {
        return std::ranges::any_of(ACandidates, [&](const Point3d &c) { return Vector3d(AP, c).norm() <= ATol; });
    }
} // namespace

TEST_CASE("quad_block_classification_snaps_corners_edges_and_face_onto_geom_model", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    // Corners deliberately offset from the real square corners, but within tolerance.
    const std::array<Point3d, 4> corners = {
        Point3d(0.01, 0.0, 0.0), Point3d(0.99, 0.02, 0.0), Point3d(1.0, 0.98, 0.0), Point3d(0.02, 1.0, 0.0)};
    auto face = blocking.create_quad_block(corners);

    blocking.classify(0.05);

    const std::vector<Point3d> expected_vertex_positions = {
        Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)};

    std::vector<Int> found_vertex_tags;
    for (auto it = blocking.cmap().attributes<0>().begin(), itend = blocking.cmap().attributes<0>().end(); it != itend;
         ++it) {
        REQUIRE(it->info().geom_targets.size() == 1);
        REQUIRE(it->info().geom_targets[0].first == GroupDim::Dim0);
        const Int tag = it->info().geom_targets[0].second;
        found_vertex_tags.push_back(tag);

        const auto *vertex = geom.vertex_by_tag(tag);
        REQUIRE(vertex != nullptr);
        REQUIRE(it->info().point == vertex->closest_point(Point3d(0, 0, 0)));
    }
    std::ranges::sort(found_vertex_tags);
    REQUIRE(found_vertex_tags == std::vector<Int>{1, 2, 3, 4});

    std::vector<Int> found_curve_tags;
    for (auto it = blocking.cmap().attributes<1>().begin(), itend = blocking.cmap().attributes<1>().end(); it != itend;
         ++it) {
        REQUIRE(it->info().geom_targets.size() == 1);
        REQUIRE(it->info().geom_targets[0].first == GroupDim::Dim1);
        found_curve_tags.push_back(it->info().geom_targets[0].second);
    }
    std::ranges::sort(found_curve_tags);
    REQUIRE(found_curve_tags == std::vector<Int>{10, 11, 12, 13});

    REQUIRE(face->info().geom_targets.size() == 1);
    REQUIRE(face->info().geom_targets[0] == std::pair{GroupDim::Dim2, Int(20)});

    // The rebuilt surface's 4 corners must still be (a permutation of) the 4 real vertex positions.
    REQUIRE(is_near_any(face->info().surface.value(0.0, 0.0), expected_vertex_positions, 1e-9));
    REQUIRE(is_near_any(face->info().surface.value(1.0, 0.0), expected_vertex_positions, 1e-9));
    REQUIRE(is_near_any(face->info().surface.value(0.0, 1.0), expected_vertex_positions, 1e-9));
    REQUIRE(is_near_any(face->info().surface.value(1.0, 1.0), expected_vertex_positions, 1e-9));
}

TEST_CASE("cubic_quad_block_edge_classified_onto_bent_curve_bulges_and_propagates_to_face", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model_with_bent_bottom();
    Blocking<FacetedGeometry> blocking(geom, 3);

    const std::array<Point3d, 4> corners = {
        Point3d(0.01, 0.0, 0.0), Point3d(0.99, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)};
    auto face = blocking.create_quad_block(corners);

    // Curve/surface tolerance: the bottom edge's own pre-classification straight midpoint is
    // ~0.31 from the bent reference curve it must still be recognized as classifiable onto, while
    // staying tight enough to exclude the left/right edges' own straight midpoints, which also
    // happen to pass within ~0.39 of the bend (0.6 for the top edge, safely excluded either way).
    blocking.classify(0.05, 0.35);

    // Find the edge classified onto the bent bottom curve (tag 10).
    using Edge = Blocking<FacetedGeometry>::Edge;
    Edge bottom_edge{};
    bool found = false;
    for (auto it = blocking.cmap().attributes<1>().begin(), itend = blocking.cmap().attributes<1>().end(); it != itend;
         ++it) {
        if (!it->info().geom_targets.empty() && it->info().geom_targets[0] == std::pair{GroupDim::Dim1, Int(10)}) {
            bottom_edge = it;
            found = true;
        }
    }
    REQUIRE(found);

    // The refit edge must genuinely bulge toward the bent curve's own apex, not stay straight.
    const Point3d mid = bottom_edge->info().curve.value(0.5);
    REQUIRE(mid.y() > 0.1);

    // Its endpoints must still be exactly its own (classified) corner nodes.
    REQUIRE(bottom_edge->info().curve.control_points()[0].y() == Approx(0.0).margin(1e-9));
    REQUIRE(bottom_edge->info().curve.control_points()[3].y() == Approx(0.0).margin(1e-9));

    // The curvature must propagate identically into the face's own stored surface: by construction
    // (see CoonsPatch.h), a Coons surface's boundary row/column exactly reproduces its boundary
    // edge curve — on whichever of the face's 4 boundaries happens to be the bottom edge (the
    // face's own local (u,v) frame, walked from its own dart, isn't assumed to start at any
    // particular physical corner).
    const std::array<Point3d, 4> face_boundary_mids = {face->info().surface.value(0.5, 0.0),
                                                       face->info().surface.value(0.5, 1.0),
                                                       face->info().surface.value(0.0, 0.5),
                                                       face->info().surface.value(1.0, 0.5)};
    REQUIRE(is_near_any(mid, {face_boundary_mids.begin(), face_boundary_mids.end()}, 1e-9));
}

TEST_CASE("node_near_untagged_curve_junction_classifies_onto_the_group_of_both_curves", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_curve_junction_geom_model();
    REQUIRE(geom.nb_vertices() == 0);
    Blocking<FacetedGeometry> blocking(geom);

    // A small quad whose corner 0 sits just off the (untagged) junction of curves 100 and 101.
    const std::array<Point3d, 4> corners = {
        Point3d(0.01, 0.01, 0.0), Point3d(0.5, 0.0, 0.0), Point3d(0.5, 0.5, 0.0), Point3d(0.0, 0.5, 0.0)};
    blocking.create_quad_block(corners);

    blocking.classify(0.02, 0.05);

    // The junction node is the one with 2 classification targets (no vertex exists to snap to, so
    // its search falls through to dimension 1 and finds both curves) — its position moves off the
    // raw (0.01,0.01,0) corner once classify() projects it onto the nearer of the 2 curves.
    bool found_junction_node = false;
    for (auto it = blocking.cmap().attributes<0>().begin(), itend = blocking.cmap().attributes<0>().end(); it != itend;
         ++it) {
        if (it->info().geom_targets.size() == 2) {
            found_junction_node = true;
            std::vector<Int> tags{it->info().geom_targets[0].second, it->info().geom_targets[1].second};
            std::ranges::sort(tags);
            REQUIRE(tags == std::vector<Int>{100, 101});
            for (const auto &target : it->info().geom_targets) {
                REQUIRE(target.first == GroupDim::Dim1);
            }
        }
    }
    REQUIRE(found_junction_node);
}

TEST_CASE("hex_block_geometry_survives_classify_rebuild_when_nothing_matches", "[BlockTestSuite]") {
    // Regression test for classify_and_rebuild_block()'s edge/face-role re-derivation: unlike
    // create_hex_block() (which builds its 6 face surfaces directly, in one pass, from a known
    // corner array), classify() must re-derive the same roles purely combinatorially, from
    // whichever dart the Block attribute happens to report as its own — not necessarily the same
    // dart (up to the hex's own rotation/reflection symmetry) used at creation, so the rebuilt
    // volume's 8 corners may come out as a *relabeled* permutation of the original 8, not
    // necessarily at the same (u,v,w) as before. What must not change is the actual shape: the set
    // of 8 corners and the centroid.
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    const std::array<Point3d, 8> corners = {Point3d(0.0, 0.0, 0.0),
                                            Point3d(1.0, 0.0, 0.0),
                                            Point3d(1.0, 1.0, 0.0),
                                            Point3d(0.0, 1.0, 0.0),
                                            Point3d(0.0, 0.0, 1.0),
                                            Point3d(1.0, 0.0, 1.0),
                                            Point3d(1.0, 1.0, 1.0),
                                            Point3d(0.0, 1.0, 1.0)};
    auto block = blocking.create_hex_block(corners);

    // The square fixture's geometry lives near the origin in the z=0 plane; use a minuscule
    // tolerance so nothing classifies (classify()'s rebuild path still runs unconditionally).
    blocking.classify(1e-12);

    const std::vector<Point3d> corner_set(corners.begin(), corners.end());
    for (double u : {0.0, 1.0}) {
        for (double v : {0.0, 1.0}) {
            for (double w : {0.0, 1.0}) {
                REQUIRE(is_near_any(block->info().volume.value(u, v, w), corner_set, 1e-9));
            }
        }
    }
    const Point3d center = block->info().volume.value(0.5, 0.5, 0.5);
    REQUIRE(center.x() == Approx(0.5));
    REQUIRE(center.y() == Approx(0.5));
    REQUIRE(center.z() == Approx(0.5));
}

TEST_CASE("hex_block_classification_against_real_multi_volume_gmsh_file", "[BlockTestSuite]") {
    std::string dir(TEST_SAMPLES_DIR);
    const auto path = dir + "/two_cubes.msh";
    const FacetedGeometry geom(path);
    Blocking<FacetedGeometry> blocking(geom);

    // Anchor one corner exactly at a real model vertex (queried, not guessed) plus a tiny offset.
    const Point3d real_vertex_pos = geom.vertices()[0].closest_point(Point3d(0, 0, 0));
    const Int expected_tag = geom.vertices()[0].entity_tag();
    const Vector3d off(0.001, 0.001, 0.001);

    const std::array<Point3d, 8> corners = {real_vertex_pos + off,
                                            real_vertex_pos + off + Vector3d(0.3, 0, 0),
                                            real_vertex_pos + off + Vector3d(0.3, 0.3, 0),
                                            real_vertex_pos + off + Vector3d(0, 0.3, 0),
                                            real_vertex_pos + off + Vector3d(0, 0, 0.3),
                                            real_vertex_pos + off + Vector3d(0.3, 0, 0.3),
                                            real_vertex_pos + off + Vector3d(0.3, 0.3, 0.3),
                                            real_vertex_pos + off + Vector3d(0, 0.3, 0.3)};
    auto block = blocking.create_hex_block(corners);

    blocking.classify(0.01);
    REQUIRE(blocking.is_valid_topology());

    bool found_snapped_vertex = false;
    for (auto it = blocking.cmap().attributes<0>().begin(), itend = blocking.cmap().attributes<0>().end(); it != itend;
         ++it) {
        if (!it->info().geom_targets.empty() && it->info().geom_targets[0].first == GroupDim::Dim0 &&
            it->info().geom_targets[0].second == expected_tag) {
            found_snapped_vertex = true;
            REQUIRE(it->info().point == real_vertex_pos);
        }
    }
    REQUIRE(found_snapped_vertex);

    // The block itself searches only dimension 3 (volumes); FacetedVolume's point-membership is a
    // documented stub (always distance 0), so every volume in the file is a valid target.
    REQUIRE(block->info().geom_targets.size() == geom.nb_volumes());
    for (const auto &target : block->info().geom_targets) {
        REQUIRE(target.first == GroupDim::Dim3);
    }

    auto mesh = blocking.to_mesh(10);
    REQUIRE(mesh.nb_cells() == 1000);
}
