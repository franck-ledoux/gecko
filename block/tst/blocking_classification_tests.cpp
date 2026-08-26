#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>
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
     * @brief A faceted unit sphere centred on the origin: one tagged surface (tag 1) of 2048
     * triangles, obtained by subdividing an octahedron 4 times and pushing every node out to radius
     * 1. No tagged vertices or curves — nothing here should ever classify below dimension 2.
     *
     * Doubly curved on purpose. A cylinder is developable along its axis, so a Coons patch spanning
     * one already follows it closely and cannot tell a face that was fitted to the geometry from one
     * that merely interpolates its 4 boundary curves. A sphere separates the 2 by a factor of 4.
     */
    FacetedGeometry make_sphere_geom_model() {
        std::array<Point3d, 6> axis = {Point3d(1, 0, 0),
                                       Point3d(-1, 0, 0),
                                       Point3d(0, 1, 0),
                                       Point3d(0, -1, 0),
                                       Point3d(0, 0, 1),
                                       Point3d(0, 0, -1)};
        std::vector<std::array<Point3d, 3>> tris;
        for (int a = 0; a < 2; ++a) {
            for (int b = 2; b < 4; ++b) {
                for (int c = 4; c < 6; ++c) {
                    tris.push_back({axis[static_cast<std::size_t>(a)],
                                    axis[static_cast<std::size_t>(b)],
                                    axis[static_cast<std::size_t>(c)]});
                }
            }
        }
        const auto on_sphere = [](const Point3d &AP) {
            const double n = std::sqrt(AP.x() * AP.x() + AP.y() * AP.y() + AP.z() * AP.z());
            return Point3d(AP.x() / n, AP.y() / n, AP.z() / n);
        };
        for (int pass = 0; pass < 4; ++pass) {
            std::vector<std::array<Point3d, 3>> finer;
            finer.reserve(tris.size() * 4);
            for (const auto &t : tris) {
                const Point3d m01 = on_sphere(t[0] + Vector3d(t[0], t[1]) * 0.5);
                const Point3d m12 = on_sphere(t[1] + Vector3d(t[1], t[2]) * 0.5);
                const Point3d m20 = on_sphere(t[2] + Vector3d(t[2], t[0]) * 0.5);
                finer.push_back({t[0], m01, m20});
                finer.push_back({m01, t[1], m12});
                finer.push_back({m20, m12, t[2]});
                finer.push_back({m01, m12, m20});
            }
            tris.swap(finer);
        }

        SimplicialMesh mesh;
        GroupRegistry groups;
        auto surf_group = groups.add_group("Sphere", GroupDim::Dim2);
        auto &face_group = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &face_entity = mesh.add_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
        std::map<std::array<long, 3>, NodeId> seen;
        const auto node_of = [&](const Point3d &AP) {
            // Keyed on rounded coordinates so the 3 triangles meeting at a node share it: the
            // subdivision recomputes each midpoint once per incident triangle.
            const std::array<long, 3> key{
                std::lround(AP.x() * 1e9), std::lround(AP.y() * 1e9), std::lround(AP.z() * 1e9)};
            const auto it = seen.find(key);
            if (it != seen.end()) return it->second;
            const NodeId id = mesh.add_node(AP.x(), AP.y(), AP.z());
            seen.emplace(key, id);
            return id;
        };
        for (const auto &t : tris) {
            const auto f = mesh.add_face(node_of(t[0]), node_of(t[1]), node_of(t[2]));
            face_group[f.value] = surf_group;
            face_entity[f.value] = 1;
        }

        const auto path = (std::filesystem::temp_directory_path() / "gecko_sphere_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }

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

TEST_CASE("a_face_classified_on_a_curved_surface_has_its_interior_on_it_too", "[BlockTestSuite]") {
    // A face's surface is a Coons patch of its 4 boundary curves, and a Coons patch of 4 curves that
    // lie on a sphere does not lie on the sphere — it interpolates the boundary and then blends
    // straight across the middle. So classifying a face onto a surface has to *pull* its interior
    // onto that surface, the way classifying an edge onto a curve pulls its interior onto the curve.
    // Without that the edges of a block sit on the model while the middle of every face hangs a
    // fifth of the radius inside it, and raising the order does not help, because nothing in the
    // construction refers to the model at all.
    const FacetedGeometry geom = make_sphere_geom_model();

    const double h = 1.0 / std::sqrt(3.0); // the cube inscribed in the unit sphere
    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_hex_block({Point3d(-h, -h, -h),
                               Point3d(h, -h, -h),
                               Point3d(h, h, -h),
                               Point3d(-h, h, -h),
                               Point3d(-h, -h, h),
                               Point3d(h, -h, h),
                               Point3d(h, h, h),
                               Point3d(-h, h, h)});
    // Nothing is tagged below dimension 2, so only the surface tolerance matters; it is opened wide
    // enough to reach the sphere from the inscribed cube's corners.
    blocking.classify(1e-9, 1e-9, 2.0);

    // How far a point is from the sphere. Measured against the analytic radius rather than the
    // faceted projection, so the fixture's own faceting shows up as a small constant on every row
    // instead of being silently absorbed.
    const auto off_sphere = [](const Point3d &AP) {
        return std::abs(std::sqrt(AP.x() * AP.x() + AP.y() * AP.y() + AP.z() * AP.z()) - 1.0);
    };

    auto &map = blocking.cmap();
    double worst_corner = 0.0;
    for (auto it = map.attributes<0>().begin(), end = map.attributes<0>().end(); it != end; ++it) {
        worst_corner = std::max(worst_corner, off_sphere(it->info().point));
    }
    double worst_edge = 0.0;
    for (auto it = map.attributes<1>().begin(), end = map.attributes<1>().end(); it != end; ++it) {
        for (int i = 1; i < 8; ++i) {
            worst_edge = std::max(worst_edge, off_sphere(it->info().curve.value(i / 8.0)));
        }
    }
    int checked = 0;
    double worst_face = 0.0;
    for (auto it = map.attributes<2>().begin(), end = map.attributes<2>().end(); it != end; ++it) {
        const auto &targets = it->info().geom_targets;
        REQUIRE(!targets.empty());
        REQUIRE(targets.front().first == GroupDim::Dim2);
        ++checked;
        for (int i = 1; i < 8; ++i) {
            for (int j = 1; j < 8; ++j) {
                worst_face = std::max(worst_face, off_sphere(it->info().surface.value(i / 8.0, j / 8.0)));
            }
        }
    }
    REQUIRE(checked == 6);

    // The corners are projected onto the geometry outright, so they are off the analytic sphere by
    // the fixture's faceting alone. The edges are fitted through points taken on it and carry, on
    // top of that, the error of a cubic through 4 points of a 70-degree arc — which a cube's edge on
    // a sphere is.
    REQUIRE(worst_corner < 5e-3);
    REQUIRE(worst_edge < 2e-2);

    // The interior lands within about twice the edges' own error. A plain Coons patch of the very
    // same 4 boundary curves misses by 0.18 here — more than a sixth of the radius — so this bound
    // sits between the 2, wide enough not to be brittle and tight enough that dropping the refit
    // fails it.
    REQUIRE(worst_face < 9e-2);

    // And the block's own volume follows its faces rather than re-deriving a Coons patch of its
    // edges, which would put its boundary where its faces are not.
    for (auto it = map.attributes<3>().begin(), end = map.attributes<3>().end(); it != end; ++it) {
        for (int i = 0; i <= 4; ++i) {
            for (int j = 0; j <= 4; ++j) {
                REQUIRE(off_sphere(it->info().volume.value(i / 4.0, j / 4.0, 0.0)) < 9e-2);
                REQUIRE(off_sphere(it->info().volume.value(i / 4.0, j / 4.0, 1.0)) < 9e-2);
                REQUIRE(off_sphere(it->info().volume.value(i / 4.0, 0.0, j / 4.0)) < 9e-2);
                REQUIRE(off_sphere(it->info().volume.value(i / 4.0, 1.0, j / 4.0)) < 9e-2);
                REQUIRE(off_sphere(it->info().volume.value(0.0, i / 4.0, j / 4.0)) < 9e-2);
                REQUIRE(off_sphere(it->info().volume.value(1.0, i / 4.0, j / 4.0)) < 9e-2);
            }
        }
    }
}

TEST_CASE("collapsing_a_sheet_keeps_the_more_constrained_of_the_2_corners_it_merges", "[BlockTestSuite]") {
    // Collapsing a sheet merges 2 corners into 1, and one of the 2 classifications has to go. Here
    // the 2 sides are as different as they get: one corner sits on a model *vertex*, the other only
    // on a boundary *curve*. The vertex is the more constrained of the 2 and lies on that curve, so
    // it is the one that survives — and the merged corner lands exactly on the square's corner
    // rather than drifting somewhere along its bottom edge.
    //
    // Deciding this by traversal order instead — which is what merging attributes does if left to
    // itself — would pull the block structure off the model's own corner, silently.
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(0.5, 0.0, 0.0), Point3d(0.5, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
    blocking.create_quad_block(
        {Point3d(0.5, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.5, 1.0, 0.0)});
    blocking.build_connectivity();
    blocking.classify(1e-6);

    auto &map = blocking.cmap();
    const auto corner_at = [&](const Point3d &AP) {
        for (auto it = map.attributes<0>().begin(), end = map.attributes<0>().end(); it != end; ++it) {
            if (Vector3d(it->info().point, AP).norm() < 1e-9) return it;
        }
        return map.attributes<0>().end();
    };
    const auto dim_at = [&](const Point3d &AP) {
        const auto n = corner_at(AP);
        REQUIRE(n != map.attributes<0>().end());
        REQUIRE(!n->info().geom_targets.empty());
        return n->info().geom_targets.front().first;
    };

    // The starting point: the square's own corner is on a model vertex, the midpoint of its bottom
    // edge only on the bottom curve.
    REQUIRE(dim_at(Point3d(0.0, 0.0, 0.0)) == GroupDim::Dim0);
    REQUIRE(dim_at(Point3d(0.5, 0.0, 0.0)) == GroupDim::Dim1);

    // The sheet through the left quad's bottom edge is that quad's 2 horizontal edges; collapsing it
    // merges (0,0,0) with (0.5,0,0) and (0,1,0) with (0.5,1,0).
    const auto target = [&]() {
        for (auto it = map.attributes<1>().begin(), end = map.attributes<1>().end(); it != end; ++it) {
            const auto d = it->dart();
            const Point3d a = map.attribute<0>(d)->info().point;
            const Point3d b = map.attribute<0>(map.beta<1>(d))->info().point;
            if (Vector3d(a, Point3d(0.0, 0.0, 0.0)).norm() < 1e-9 &&
                Vector3d(b, Point3d(0.5, 0.0, 0.0)).norm() < 1e-9) {
                return it;
            }
            if (Vector3d(b, Point3d(0.0, 0.0, 0.0)).norm() < 1e-9 &&
                Vector3d(a, Point3d(0.5, 0.0, 0.0)).norm() < 1e-9) {
                return it;
            }
        }
        FAIL("no bottom-left edge");
        return map.attributes<1>().end();
    }();
    REQUIRE(blocking.delete_sheet(target, 1e-6));
    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.nb_cells<2>() == 1);

    // The merged corner went to the model vertex, not to the middle of the 2 and not along the
    // curve: the more constrained side won, and it kept its own classification.
    REQUIRE(corner_at(Point3d(0.0, 0.0, 0.0)) != map.attributes<0>().end());
    REQUIRE(corner_at(Point3d(0.0, 1.0, 0.0)) != map.attributes<0>().end());
    REQUIRE(dim_at(Point3d(0.0, 0.0, 0.0)) == GroupDim::Dim0);
    REQUIRE(dim_at(Point3d(0.0, 1.0, 0.0)) == GroupDim::Dim0);
    REQUIRE(corner_at(Point3d(0.5, 0.0, 0.0)) == map.attributes<0>().end());

    // And the surviving quad is the whole square again, on its 4 model vertices.
    REQUIRE(blocking.nb_cells<0>() == 4);
    REQUIRE(dim_at(Point3d(1.0, 0.0, 0.0)) == GroupDim::Dim0);
    REQUIRE(dim_at(Point3d(1.0, 1.0, 0.0)) == GroupDim::Dim0);
}

TEST_CASE("collapsing_is_refused_when_it_would_merge_2_different_model_vertices", "[BlockTestSuite]") {
    // The one thing a collapse cannot do. The square's bottom edge joins 2 corners sitting on 2
    // *different* model vertices, and merging them would leave one of those vertices with no corner
    // of the block structure on it — nothing else in the blocking records where it was, and no later
    // classify() puts one back. Refused, and nothing moves.
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
    blocking.build_connectivity();
    blocking.classify(1e-6);

    auto &map = blocking.cmap();
    // All 4 corners are on the square's own vertices, so every edge of it joins 2 different ones.
    for (auto it = map.attributes<0>().begin(), end = map.attributes<0>().end(); it != end; ++it) {
        REQUIRE(!it->info().geom_targets.empty());
        REQUIRE(it->info().geom_targets.front().first == GroupDim::Dim0);
    }

    for (auto it = map.attributes<1>().begin(), end = map.attributes<1>().end(); it != end; ++it) {
        REQUIRE_FALSE(blocking.delete_sheet(it, 1e-6));
    }

    // Untouched: not a corner moved, not a cell gone.
    REQUIRE(blocking.nb_cells<2>() == 1);
    REQUIRE(blocking.nb_cells<1>() == 4);
    REQUIRE(blocking.nb_cells<0>() == 4);
    for (const Point3d &corner :
         {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)}) {
        bool found = false;
        for (auto it = map.attributes<0>().begin(), end = map.attributes<0>().end(); it != end; ++it) {
            if (Vector3d(it->info().point, corner).norm() < 1e-12) found = true;
        }
        REQUIRE(found);
    }

    // Cut it in 2 and the *inner* sheet becomes collapsible again: its edges join a corner on a
    // vertex to one merely on a curve, and the vertex survives that. Only the vertex-to-vertex
    // pairing is the problem, not classification as such.
    const auto vertical = [&]() {
        for (auto it = map.attributes<1>().begin(), end = map.attributes<1>().end(); it != end; ++it) {
            const auto d = it->dart();
            const Point3d a = map.attribute<0>(d)->info().point;
            const Point3d b = map.attribute<0>(map.beta<1>(d))->info().point;
            if (std::abs(a.x() - b.x()) < 1e-12) return it;
        }
        FAIL("no vertical edge");
        return map.attributes<1>().end();
    }();
    (void)vertical;
    const auto horizontal = [&]() {
        for (auto it = map.attributes<1>().begin(), end = map.attributes<1>().end(); it != end; ++it) {
            const auto d = it->dart();
            const Point3d a = map.attribute<0>(d)->info().point;
            const Point3d b = map.attribute<0>(map.beta<1>(d))->info().point;
            if (std::abs(a.y() - b.y()) < 1e-12) return it;
        }
        FAIL("no horizontal edge");
        return map.attributes<1>().end();
    }();
    REQUIRE(blocking.cut_sheet(horizontal, 0.5));
    REQUIRE(blocking.nb_cells<2>() == 2);

    bool collapsed = false;
    for (auto it = map.attributes<1>().begin(), end = map.attributes<1>().end(); it != end; ++it) {
        if (blocking.delete_sheet(it, 1e-6)) {
            collapsed = true;
            break;
        }
    }
    REQUIRE(collapsed);
    REQUIRE(blocking.nb_cells<2>() == 1);
    REQUIRE(blocking.is_valid_topology());
}

TEST_CASE("cutting_a_standalone_quad_block_splits_it_whichever_way_its_sheet_runs", "[BlockTestSuite]") {
    // A boundary edge of a standalone quad block has exactly one dart — nothing is 2-sewn to it — so
    // only one of its 2 endpoints is the source of any dart of it at all. The cut used to ask for a
    // dart starting at the end its sheet measures from, and on a 2D block that ask has no answer
    // half the time: whether it aborted came down to which endpoint the sheet happened to start at.
    //
    // A hex hides this completely. Its edges carry darts in several faces, in both directions, so
    // the ask always succeeds — which is why every cut test until now passed.
    const FacetedGeometry geom = make_square_geom_model();

    for (const double param : {0.5, 0.25, 0.75}) {
        Blocking<FacetedGeometry> blocking(geom, 3);
        blocking.create_quad_block(
            {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
        blocking.build_connectivity();
        blocking.classify(1e-6);

        // Every edge in turn, so neither orientation of a sheet goes untried.
        auto &map = blocking.cmap();
        int cuts = 0;
        for (auto it = map.attributes<1>().begin(), end = map.attributes<1>().end(); it != end; ++it) {
            Blocking<FacetedGeometry> one(geom, 3);
            one.create_quad_block(
                {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
            one.build_connectivity();
            one.classify(1e-6);

            auto &m = one.cmap();
            auto target = m.attributes<1>().begin();
            std::advance(target, cuts);
            REQUIRE(one.cut_sheet(target, param));
            REQUIRE(one.is_valid_topology());
            REQUIRE(one.nb_cells<2>() == 2);
            REQUIRE(one.nb_cells<1>() == 7);
            REQUIRE(one.nb_cells<0>() == 6);

            // The cut lands where it was asked to, on the right side of the edge it was aimed at.
            const auto mesh = one.to_mesh(1);
            double lo = 1e30;
            double hi = -1e30;
            for (UInt i = 0; i < mesh.nb_nodes(); ++i) {
                const Point3d &q = mesh.node(NodeId{i});
                lo = std::min({lo, q.x(), q.y()});
                hi = std::max({hi, q.x(), q.y()});
            }
            REQUIRE(lo == Approx(0.0).margin(1e-12));
            REQUIRE(hi == Approx(1.0).margin(1e-12));
            ++cuts;
        }
        REQUIRE(cuts == 4);
    }
}

TEST_CASE("a_classified_face_follows_its_surface_as_closely_as_its_own_edges_do", "[BlockTestSuite]") {
    // Pulling a face's interior onto its surface is not enough on its own: pinning it at the
    // (degree+1)^2 parameters of an interpolation says nothing about what happens *between* them,
    // and between them is where the error lives. Measured on the sphere at degree 3, the interpolant
    // sat 4e-03 off at its own parameters — the fixture's own faceting, so exact — and 2.2e-02 off
    // halfway between them.
    //
    // Fitted over far more samples than there are unknowns instead, the same way an edge on a curve
    // is fitted, a face reaches what its own edges reach. Stated at degree 5, because at degree 3 the
    // limit is the bicubic itself: a patch of that degree cannot follow a 70-degree spherical cap
    // more closely than about 5% of the radius, whichever way its interior is chosen.
    const FacetedGeometry geom = make_sphere_geom_model();

    const double h = 1.0 / std::sqrt(3.0);
    Blocking<FacetedGeometry> blocking(geom, 5);
    blocking.create_hex_block({Point3d(-h, -h, -h),
                               Point3d(h, -h, -h),
                               Point3d(h, h, -h),
                               Point3d(-h, h, -h),
                               Point3d(-h, -h, h),
                               Point3d(h, -h, h),
                               Point3d(h, h, h),
                               Point3d(-h, h, h)});
    blocking.classify(1e-9, 1e-9, 2.0);

    const auto off_sphere = [](const Point3d &AP) {
        return std::abs(std::sqrt(AP.x() * AP.x() + AP.y() * AP.y() + AP.z() * AP.z()) - 1.0);
    };

    auto &map = blocking.cmap();
    double worst_edge = 0.0;
    for (auto it = map.attributes<1>().begin(), end = map.attributes<1>().end(); it != end; ++it) {
        for (int i = 1; i < 16; ++i) {
            worst_edge = std::max(worst_edge, off_sphere(it->info().curve.value(i / 16.0)));
        }
    }

    double worst_face = 0.0;
    int checked = 0;
    for (auto it = map.attributes<2>().begin(), end = map.attributes<2>().end(); it != end; ++it) {
        REQUIRE(!it->info().geom_targets.empty());
        REQUIRE(it->info().geom_targets.front().first == GroupDim::Dim2);
        ++checked;
        for (int i = 1; i < 16; ++i) {
            for (int j = 1; j < 16; ++j) {
                worst_face = std::max(worst_face, off_sphere(it->info().surface.value(i / 16.0, j / 16.0)));
            }
        }
    }
    REQUIRE(checked == 6);

    // The floor both are working against is the fixture's own faceting, which the corners sit at.
    double worst_corner = 0.0;
    for (auto it = map.attributes<0>().begin(), end = map.attributes<0>().end(); it != end; ++it) {
        worst_corner = std::max(worst_corner, off_sphere(it->info().point));
    }
    REQUIRE(worst_corner < 5e-3);
    REQUIRE(worst_edge < 1e-2);

    // The point of the whole thing: a face is no further off its surface than the edges bounding it.
    // Interpolating instead left it at 1.05e-02 here, half again as far.
    REQUIRE(worst_face < 1.5 * worst_edge);
}

TEST_CASE("a_curve_classified_edge_follows_its_curve_at_every_degree", "[BlockTestSuite]") {
    // The regression case for issue #48: an edge classified on a model curve has to have its own
    // sampled points sit close to that curve, not only its 2 endpoints. Measured across every degree
    // rather than at just one, because the defect this covers was invisible at degree 3 alone —
    // degree 2 collapsed a shared array slot onto itself and threw away half a least-squares solve,
    // and degree 4 and up quietly pinned their extra control points to the straight chord instead of
    // fitting them, both of which left the worst-followed degree *further* from the curve than the
    // straight-edge (degree 1) case it was supposed to improve on.
    std::string dir(TEST_SAMPLES_DIR);
    const FacetedGeometry geom(dir + "/cylinder.msh");

    double lo[3] = {1e30, 1e30, 1e30};
    double hi[3] = {-1e30, -1e30, -1e30};
    for (UInt i = 0; i < geom.mesh().nb_nodes(); ++i) {
        const Point3d &p = geom.mesh().node(NodeId{i});
        const std::array<double, 3> c{p.x(), p.y(), p.z()};
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], c[k]);
            hi[k] = std::max(hi[k], c[k]);
        }
    }
    const double extent = std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});

    double previous = 0.0;
    for (const std::size_t degree : {std::size_t{2}, std::size_t{3}, std::size_t{4}, std::size_t{5}, std::size_t{7}}) {
        Blocking<FacetedGeometry> blocking(geom, degree);
        blocking.create_hex_block({Point3d(lo[0], lo[1], lo[2]),
                                   Point3d(hi[0], lo[1], lo[2]),
                                   Point3d(hi[0], hi[1], lo[2]),
                                   Point3d(lo[0], hi[1], lo[2]),
                                   Point3d(lo[0], lo[1], hi[2]),
                                   Point3d(hi[0], lo[1], hi[2]),
                                   Point3d(hi[0], hi[1], hi[2]),
                                   Point3d(lo[0], hi[1], hi[2])});
        blocking.classify(0.3);

        double worst = 0.0;
        int checked = 0;
        auto &map = blocking.cmap();
        for (auto it = map.attributes<1>().begin(), end = map.attributes<1>().end(); it != end; ++it) {
            if (it->info().geom_targets.empty() || it->info().geom_targets.front().first != GroupDim::Dim1) continue;
            const auto *curve = geom.curve_by_tag(it->info().geom_targets.front().second);
            REQUIRE(curve != nullptr);
            ++checked;
            for (int i = 1; i < 16; ++i) {
                worst = std::max(worst, curve->distance(it->info().curve.value(i / 16.0)));
            }
        }
        REQUIRE(checked > 0);

        // Loose in absolute terms — this is a coarse faceted cylinder, not a smooth analytic one —
        // but tight enough to catch either defect: both used to push a curve *off* the model by
        // several percent of its own extent, well past what any of these degrees now reaches.
        REQUIRE(worst < 0.01 * extent);
        INFO("degree " << degree << " worst deviation " << worst << " (" << 100.0 * worst / extent << "% of extent)");
        // Not strictly monotone — a coarse facet can favour one degree over its neighbour by chance
        // — but a fit that keeps adding control points should not be getting *worse* by an order of
        // magnitude, which is what the bug did between degree 2 and degree 3 before this fix.
        if (previous > 0.0) REQUIRE(worst < 5.0 * previous);
        previous = worst;
    }
}

TEST_CASE("an_edge_on_a_surface_lies_in_one_plane_and_not_just_somewhere_on_it", "[BlockTestSuite]") {
    // The half of #48 about edges. A surface holds infinitely many curves between one pair of its
    // points, so "classified on a surface" says where the edge lands but not which way it goes —
    // and projecting each sample separately traces whichever curve the projection happens to give.
    //
    // The edge is pinned instead to the surface's section by the plane through its 2 ends containing
    // the surface's normal. What that buys is measured here: the edge lies in that plane, rather
    // than wandering a quarter of a percent of the model out of it as it did before.
    std::string dir(TEST_SAMPLES_DIR);
    const FacetedGeometry geom(dir + "/cylinder.msh");

    double lo[3] = {1e30, 1e30, 1e30};
    double hi[3] = {-1e30, -1e30, -1e30};
    for (UInt i = 0; i < geom.mesh().nb_nodes(); ++i) {
        const Point3d &p = geom.mesh().node(NodeId{i});
        const std::array<double, 3> c{p.x(), p.y(), p.z()};
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], c[k]);
            hi[k] = std::max(hi[k], c[k]);
        }
    }

    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_hex_block({Point3d(lo[0], lo[1], lo[2]),
                               Point3d(hi[0], lo[1], lo[2]),
                               Point3d(hi[0], hi[1], lo[2]),
                               Point3d(lo[0], hi[1], lo[2]),
                               Point3d(lo[0], lo[1], hi[2]),
                               Point3d(hi[0], lo[1], hi[2]),
                               Point3d(hi[0], hi[1], hi[2]),
                               Point3d(lo[0], hi[1], hi[2])});
    blocking.classify(0.3);

    auto &map = blocking.cmap();
    int checked = 0;
    for (auto it = map.attributes<1>().begin(), end = map.attributes<1>().end(); it != end; ++it) {
        const auto &targets = it->info().geom_targets;
        if (targets.empty() || targets.front().first != GroupDim::Dim2) continue;
        const auto *surface = geom.surface_by_tag(targets.front().second);
        REQUIRE(surface != nullptr);
        ++checked;

        const auto &curve = it->info().curve;
        const Point3d p0 = curve.value(0.0);
        const Point3d p1 = curve.value(1.0);
        const Vector3d chord(p0, p1);

        // The very plane the fit is pinned to: through both ends, along the surface's own normal.
        Point3d middle = p0 + chord * 0.5;
        surface->project(middle);
        const Vector3d plane_normal = chord.cross(surface->normal(middle)).normalized();

        double off_plane = 0.0;
        double off_surface = 0.0;
        for (int i = 0; i <= 24; ++i) {
            const Point3d at = curve.value(i / 24.0);
            off_plane = std::max(off_plane, std::abs(Vector3d(p0, at).dot(plane_normal)));
            off_surface = std::max(off_surface, surface->distance(at));
        }

        // In the plane, to rounding. Projecting sample by sample left it 5.7e-03 out of it — a
        // quarter of a percent of a model 2 units across.
        REQUIRE(off_plane < 1e-6);
        // And still on the surface: the fixture's own faceting is the floor, and this stays at it.
        REQUIRE(off_surface < 1e-2);
    }
    REQUIRE(checked == 4);
}

TEST_CASE("an_edge_along_a_ruling_of_a_surface_stays_straight", "[BlockTestSuite]") {
    // The trap the section fit fell into twice, and the reason it is fitted the way it is. These 4
    // edges run along the cylinder's axis, where the surface is ruled: the right answer is a
    // straight line, so any bow at all is the fit's own invention.
    //
    // Constraining the end tangents bowed them by 5e-02 — solving for 2 tangent magnitudes on data
    // that is nearly straight is nearly singular, and the control points came out oscillating.
    // Interpolating the section at uniform parameters bowed them by 2e-02 at degree 6, which is
    // Runge's phenomenon. Least squares with the 2 ends pinned does neither.
    std::string dir(TEST_SAMPLES_DIR);
    const FacetedGeometry geom(dir + "/cylinder.msh");

    double lo[3] = {1e30, 1e30, 1e30};
    double hi[3] = {-1e30, -1e30, -1e30};
    for (UInt i = 0; i < geom.mesh().nb_nodes(); ++i) {
        const Point3d &p = geom.mesh().node(NodeId{i});
        const std::array<double, 3> c{p.x(), p.y(), p.z()};
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], c[k]);
            hi[k] = std::max(hi[k], c[k]);
        }
    }

    for (const std::size_t degree : {std::size_t{2}, std::size_t{3}, std::size_t{5}, std::size_t{7}}) {
        Blocking<FacetedGeometry> blocking(geom, degree);
        blocking.create_hex_block({Point3d(lo[0], lo[1], lo[2]),
                                   Point3d(hi[0], lo[1], lo[2]),
                                   Point3d(hi[0], hi[1], lo[2]),
                                   Point3d(lo[0], hi[1], lo[2]),
                                   Point3d(lo[0], lo[1], hi[2]),
                                   Point3d(hi[0], lo[1], hi[2]),
                                   Point3d(hi[0], hi[1], hi[2]),
                                   Point3d(lo[0], hi[1], hi[2])});
        blocking.classify(0.3);

        auto &map = blocking.cmap();
        double worst = 0.0;
        for (auto it = map.attributes<1>().begin(), end = map.attributes<1>().end(); it != end; ++it) {
            const auto &targets = it->info().geom_targets;
            if (targets.empty() || targets.front().first != GroupDim::Dim2) continue;

            const auto &curve = it->info().curve;
            const Point3d p0 = curve.value(0.0);
            const Point3d p1 = curve.value(1.0);
            const Vector3d chord(p0, p1);
            const double length = chord.norm();
            if (length <= 0.0) continue;
            const Vector3d along = chord.normalized();
            for (int i = 1; i < 24; ++i) {
                const Point3d at = curve.value(i / 24.0);
                const Vector3d from_start(p0, at);
                const Vector3d across = from_start - along * from_start.dot(along);
                worst = std::max(worst, across.norm());
            }
        }
        // A hundredth of the fixture's own faceting away from straight, at every degree. It is not
        // 0 because the surface it is cut from is a prism, not a cylinder.
        REQUIRE(worst < 1e-2);
    }
}
