#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <set>
#include <vector>

#include <gecko/block/Blocking.h>
#include <gecko/geom/FacetedGeometry.h>
#include <gecko/io/GmshMeshWriter.h>
#include <unit_test_config.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
using Catch::Approx;

using namespace gecko;

namespace {
    /** @brief A geometric model with nothing in it worth classifying onto: a single far-away
     * triangle, so that every tolerance question has the same answer — no. */
    FacetedGeometry make_far_away_geom_model() {
        SimplicialMesh mesh;
        auto n0 = mesh.add_node(100, 100, 100);
        auto n1 = mesh.add_node(101, 100, 100);
        auto n2 = mesh.add_node(100, 101, 100);

        GroupRegistry groups;
        auto surf = groups.add_group("Surf", GroupDim::Dim2);
        auto &face_group = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &face_entity = mesh.add_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
        auto f0 = mesh.add_face(n0, n1, n2);
        face_group[f0.value] = surf;
        face_entity[f0.value] = 1;

        const auto path = (std::filesystem::temp_directory_path() / "gecko_block_pillow_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }

    /** @brief A geometric model whose only features are 8 vertices, one at each corner of the unit
     * cube — so that a block fitted to it has all 8 of its own corners on a *different* model
     * vertex, which is what makes any merge between 2 of them a loss of information. */
    FacetedGeometry make_cube_vertices_geom_model() {
        SimplicialMesh mesh;
        GroupRegistry groups;
        auto vtx_group = groups.add_group("Vertices", GroupDim::Dim0);
        auto surf_group = groups.add_group("Surf", GroupDim::Dim2);

        auto &node_group = mesh.add_variable<GroupId, CellType::Node>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &node_entity = mesh.add_variable<Int, CellType::Node>(std::string(io::ENTITY_TAG_VARIABLE));
        Int tag = 1;
        for (const int i : {0, 1}) {
            for (const int j : {0, 1}) {
                for (const int k : {0, 1}) {
                    auto n = mesh.add_node(i, j, k);
                    node_group[n.value] = vtx_group;
                    node_entity[n.value] = tag++;
                }
            }
        }

        // A surface has to exist for the model to be one; it is put out of every tolerance's way.
        auto f0 = mesh.add_node(100, 100, 100);
        auto f1 = mesh.add_node(101, 100, 100);
        auto f2 = mesh.add_node(100, 101, 100);
        auto &face_group = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &face_entity = mesh.add_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
        auto f = mesh.add_face(f0, f1, f2);
        face_group[f.value] = surf_group;
        face_entity[f.value] = 1;

        const auto path = (std::filesystem::temp_directory_path() / "gecko_block_pillow_vertices.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }

    /** @brief The 8 corners of the axis-aligned box `[x0,x1] x [y0,y1] x [0,1]`, in HEX8 order. */
    std::array<Point3d, 8> box(double x0, double x1, double y0 = 0.0, double y1 = 1.0) {
        return {Point3d(x0, y0, 0),
                Point3d(x1, y0, 0),
                Point3d(x1, y1, 0),
                Point3d(x0, y1, 0),
                Point3d(x0, y0, 1),
                Point3d(x1, y0, 1),
                Point3d(x1, y1, 1),
                Point3d(x0, y1, 1)};
    }

    /** @brief The faces bounding one block. */
    std::vector<Blocking<FacetedGeometry>::Face> faces_of(Blocking<FacetedGeometry> &ABlocking,
                                                          Blocking<FacetedGeometry>::Block ABlock) {
        std::vector<Blocking<FacetedGeometry>::Face> faces;
        auto &map = ABlocking.cmap();
        for (auto it = map.one_dart_per_incident_cell<2, 3>(ABlock->dart()).begin(),
                  end = map.one_dart_per_incident_cell<2, 3>(ABlock->dart()).end();
             it != end;
             ++it) {
            faces.push_back(map.attribute<2>(it));
        }
        return faces;
    }

    /** @brief The faces of one block that have a block on their other side too. */
    std::vector<Blocking<FacetedGeometry>::Face> shared_faces_of(Blocking<FacetedGeometry> &ABlocking,
                                                                 Blocking<FacetedGeometry>::Block ABlock) {
        std::vector<Blocking<FacetedGeometry>::Face> shared;
        for (const auto f : faces_of(ABlocking, ABlock)) {
            if (!ABlocking.cmap().is_free<3>(f->dart())) shared.push_back(f);
        }
        return shared;
    }

    /** @brief What every block of the blocking measures, summed. */
    double total_volume(Blocking<FacetedGeometry> &ABlocking) {
        double total = 0.0;
        for (const double v : ABlocking.block_volumes(2)) {
            total += v;
        }
        return total;
    }

    /** @brief Whether every block came out the right way round — which is what says each one was
     * built on the right face of its own hexahedron rather than inside out. */
    bool all_positive(Blocking<FacetedGeometry> &ABlocking) {
        for (const double v : ABlocking.block_volumes(2)) {
            if (v <= 0.0) return false;
        }
        return true;
    }

    /** @brief The face whose 4 corners sit exactly at @p AWanted, in any order. */
    Blocking<FacetedGeometry>::Face face_at(Blocking<FacetedGeometry> &ABlocking,
                                            const std::vector<std::array<double, 3>> &AWanted) {
        auto &map = ABlocking.cmap();
        const auto same = [](const Point3d &AP, const std::array<double, 3> &AQ) {
            return std::abs(AP.x() - AQ[0]) < 1e-9 && std::abs(AP.y() - AQ[1]) < 1e-9 &&
                   std::abs(AP.z() - AQ[2]) < 1e-9;
        };
        for (auto it = map.attributes<2>().begin(), end = map.attributes<2>().end(); it != end; ++it) {
            int matched = 0;
            auto walk = it->dart();
            for (int c = 0; c < 4; ++c) {
                const Point3d &p = map.attribute<0>(walk)->info().point;
                for (const auto &q : AWanted) {
                    if (same(p, q)) ++matched;
                }
                walk = map.beta<1>(walk);
            }
            if (matched == 4) return it;
        }
        return nullptr;
    }

    /** @brief The corner of @p AFace sitting at @p AAt. */
    Blocking<FacetedGeometry>::Node corner_at(Blocking<FacetedGeometry> &ABlocking,
                                              Blocking<FacetedGeometry>::Face AFace,
                                              const std::array<double, 3> &AAt) {
        auto &map = ABlocking.cmap();
        auto walk = AFace->dart();
        for (int c = 0; c < 4; ++c) {
            const auto node = map.attribute<0>(walk);
            const Point3d &p = node->info().point;
            if (std::abs(p.x() - AAt[0]) < 1e-9 && std::abs(p.y() - AAt[1]) < 1e-9 && std::abs(p.z() - AAt[2]) < 1e-9) {
                return node;
            }
            walk = map.beta<1>(walk);
        }
        return nullptr;
    }

    /** @brief Whether 2 blocks share a face. */
    bool share_a_face(Blocking<FacetedGeometry> &ABlocking,
                      Blocking<FacetedGeometry>::Block AA,
                      Blocking<FacetedGeometry>::Block AB) {
        auto &map = ABlocking.cmap();
        for (auto it = map.one_dart_per_incident_cell<2, 3>(AA->dart()).begin(),
                  end = map.one_dart_per_incident_cell<2, 3>(AA->dart()).end();
             it != end;
             ++it) {
            if (map.is_free<3>(it)) continue;
            if (map.attribute<3>(map.beta<3>(it)) == AB) return true;
        }
        return false;
    }

    /** @brief A 2x2x2 grid of unit blocks, sewn. */
    void build_grid_2x2x2(Blocking<FacetedGeometry> &ABlocking) {
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    const double x = i;
                    const double y = j;
                    const double z = k;
                    ABlocking.create_hex_block({Point3d(x, y, z),
                                                Point3d(x + 1, y, z),
                                                Point3d(x + 1, y + 1, z),
                                                Point3d(x, y + 1, z),
                                                Point3d(x, y, z + 1),
                                                Point3d(x + 1, y, z + 1),
                                                Point3d(x + 1, y + 1, z + 1),
                                                Point3d(x, y + 1, z + 1)});
                }
            }
        }
        ABlocking.build_connectivity();
    }
} // namespace

TEST_CASE("pillowing_a_lone_block_wraps_it_in_a_shell_of_6", "[BlockTestSuite]") {
    // The smallest closed nappe there is: a block's own 6 faces. The block shrinks and the gap it
    // leaves is filled by 1 new block per face — 6 of them, meeting each other along the 12 edges of
    // the original, which is what makes the shell conformal rather than 6 loose slabs.
    const FacetedGeometry geom = make_far_away_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 1);
    const auto b = blocking.create_hex_block(box(0, 1));

    const auto nappe = faces_of(blocking, b);
    REQUIRE(nappe.size() == 6);
    REQUIRE(blocking.pillow(nappe, b, 0.25, 1e-9));

    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.nb_cells<3>() == 7);
    // 8 corners of the shrunk block, and the 8 of the original, which the shell's outside keeps.
    REQUIRE(blocking.nb_cells<0>() == 16);
    // 6 inner, 6 outer, and 12 between neighbouring blocks of the shell, one per edge of the cube.
    REQUIRE(blocking.nb_cells<2>() == 24);
    REQUIRE(all_positive(blocking));
    // The shell fills exactly what the shrink emptied: nothing left the original block.
    REQUIRE(total_volume(blocking) == Approx(1.0).margin(1e-9));
}

TEST_CASE("pillowing_one_shared_face_inserts_a_block_between_its_2_neighbours", "[BlockTestSuite]") {
    // A nappe of a single interior face. The side named as the inside is the only one that moves,
    // the other staying exactly where it was — so the new block sits entirely within what the shrunk
    // one gave up, and the total is unchanged.
    const FacetedGeometry geom = make_far_away_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 1);
    const auto left = blocking.create_hex_block(box(0, 1));
    blocking.create_hex_block(box(1, 2));
    blocking.build_connectivity();
    REQUIRE(blocking.nb_cells<3>() == 2);
    REQUIRE(blocking.nb_cells<0>() == 12);

    const auto nappe = shared_faces_of(blocking, left);
    REQUIRE(nappe.size() == 1);
    REQUIRE(blocking.pillow(nappe, left, 0.25, 1e-9));

    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.nb_cells<3>() == 3);
    REQUIRE(blocking.nb_cells<0>() == 16);
    REQUIRE(all_positive(blocking));
    REQUIRE(total_volume(blocking) == Approx(2.0).margin(1e-9));

    // Straight in, so the layer is exactly as thick as it was asked to be: a quarter of the edge
    // length, and no sideways drift. A shrink towards the middle of the block instead would have
    // dragged all 4 corners diagonally and sheared it.
    const auto layer = blocking.block_volumes(2);
    std::vector<double> sorted(layer.begin(), layer.end());
    std::ranges::sort(sorted);
    REQUIRE(sorted[0] == Approx(0.25).margin(1e-9));
    REQUIRE(sorted[1] == Approx(0.75).margin(1e-9));
    REQUIRE(sorted[2] == Approx(1.0).margin(1e-9));
}

TEST_CASE("pillowing_a_boundary_face_inserts_a_layer_under_it", "[BlockTestSuite]") {
    // A nappe of a single face on the boundary of the blocking. There is nothing to unsew and no
    // outside copy to take: the layer's outer corners are simply left where the boundary was, which
    // is what makes this the boundary-layer case rather than a special one.
    const FacetedGeometry geom = make_far_away_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 1);
    const auto b = blocking.create_hex_block(box(0, 1));
    const auto nappe = faces_of(blocking, b);

    REQUIRE(blocking.pillow({nappe[0]}, b, 0.25, 1e-9));
    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.nb_cells<3>() == 2);
    REQUIRE(blocking.nb_cells<0>() == 12);
    REQUIRE(all_positive(blocking));
    REQUIRE(total_volume(blocking) == Approx(1.0).margin(1e-9));
}

TEST_CASE("pillowing_a_nappe_through_the_structure_cuts_it_in_two", "[BlockTestSuite]") {
    // The other kind of nappe: not closed, but running clean through the blocking and coming out on
    // its boundary. A 3x1x1 row, cut after the first block.
    const FacetedGeometry geom = make_far_away_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 1);
    const auto first = blocking.create_hex_block(box(0, 1));
    blocking.create_hex_block(box(1, 2));
    blocking.create_hex_block(box(2, 3));
    blocking.build_connectivity();
    REQUIRE(blocking.nb_cells<3>() == 3);

    const auto nappe = shared_faces_of(blocking, first);
    REQUIRE(nappe.size() == 1);
    REQUIRE(blocking.pillow(nappe, first, 0.25, 1e-9));

    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.nb_cells<3>() == 4);
    REQUIRE(all_positive(blocking));
    REQUIRE(total_volume(blocking) == Approx(3.0).margin(1e-9));
}

TEST_CASE("pillowing_refuses_what_is_not_a_nappe", "[BlockTestSuite]") {
    // Every refusal has to leave the blocking exactly as it was, so each is checked against the
    // counts as well as against the answer.
    const FacetedGeometry geom = make_far_away_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 1);
    // A 2x2 grid, so that 2 blocks sharing a face are also joined the long way round.
    const auto b = blocking.create_hex_block(box(0, 1, 0, 1));
    blocking.create_hex_block(box(1, 2, 0, 1));
    blocking.create_hex_block(box(0, 1, 1, 2));
    blocking.create_hex_block(box(1, 2, 1, 2));
    blocking.build_connectivity();
    REQUIRE(blocking.nb_cells<3>() == 4);
    // Built up front rather than inside the section that needs it, so that every section counts the
    // same cells before and after.
    const auto quad =
        blocking.create_quad_block({Point3d(0, 0, 5), Point3d(1, 0, 5), Point3d(1, 1, 5), Point3d(0, 1, 5)});
    const auto blocks_before = blocking.nb_cells<3>();
    const auto nodes_before = blocking.nb_cells<0>();

    const auto shared = shared_faces_of(blocking, b);
    REQUIRE(shared.size() == 2);

    SECTION("a thickness outside (0,1)") {
        REQUIRE_FALSE(blocking.pillow(shared, b, 0.0, 1e-9));
        REQUIRE_FALSE(blocking.pillow(shared, b, 1.0, 1e-9));
    }
    SECTION("no faces at all") { REQUIRE_FALSE(blocking.pillow({}, b, 0.25, 1e-9)); }
    SECTION("the same face named twice") { REQUIRE_FALSE(blocking.pillow({shared[0], shared[0]}, b, 0.25, 1e-9)); }
    SECTION("a face that does not separate") {
        // One interior face of the grid, on its own. Its 2 blocks are still joined round the other
        // 2, so nothing is on either side of it and there is no gap to insert anything into.
        REQUIRE_FALSE(blocking.pillow({shared[0]}, b, 0.25, 1e-9));
    }
    SECTION("a standalone quad block, which is not a block face at all") {
        REQUIRE_FALSE(blocking.pillow({quad}, b, 0.25, 1e-9));
    }
    REQUIRE(blocking.nb_cells<3>() == blocks_before);
    REQUIRE(blocking.nb_cells<0>() == nodes_before);
}

TEST_CASE("pillowing_takes_the_inside_off_the_model_and_leaves_the_outside_on_it", "[BlockTestSuite]") {
    // What the operation does to classification, which is the half of it that is a choice. The
    // corner the nappe cuts through becomes 2: the outside one keeps everything it had, and the
    // inside one keeps only what it is still on after moving — nothing, here, having been pushed a
    // quarter of an edge into the interior.
    std::string dir(TEST_SAMPLES_DIR);
    const FacetedGeometry geom(dir + "/cylinder.msh");
    Blocking<FacetedGeometry> blocking(geom, 1);

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
    const auto b = blocking.create_hex_block({Point3d(lo[0], lo[1], lo[2]),
                                              Point3d(hi[0], lo[1], lo[2]),
                                              Point3d(hi[0], hi[1], lo[2]),
                                              Point3d(lo[0], hi[1], lo[2]),
                                              Point3d(lo[0], lo[1], hi[2]),
                                              Point3d(hi[0], lo[1], hi[2]),
                                              Point3d(hi[0], hi[1], hi[2]),
                                              Point3d(lo[0], hi[1], hi[2])});
    blocking.classify(0.3);

    // Every corner of the single block sits on a model vertex to begin with.
    auto &map = blocking.cmap();
    for (auto it = map.attributes<0>().begin(), end = map.attributes<0>().end(); it != end; ++it) {
        REQUIRE_FALSE(it->info().geom_targets.empty());
    }

    REQUIRE(blocking.pillow(faces_of(blocking, b), b, 0.25, 1e-6, 1e-3, 1e-2));
    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.nb_cells<3>() == 7);
    REQUIRE(all_positive(blocking));

    int classified = 0;
    int unclassified = 0;
    for (auto it = map.attributes<0>().begin(), end = map.attributes<0>().end(); it != end; ++it) {
        if (it->info().geom_targets.empty()) {
            ++unclassified;
        } else {
            ++classified;
        }
    }
    // The 8 that stayed are still on the model; the 8 that moved inwards are on nothing at all.
    REQUIRE(classified == 8);
    REQUIRE(unclassified == 8);
}

TEST_CASE("pillowing_an_unclassified_blocking_leaves_it_unclassified", "[BlockTestSuite]") {
    // The mirror of the case above, and the reason the rule is stated as "keeps what it is still on"
    // rather than as a proximity search: a blocking nobody classified must come out of a pillow with
    // nothing classified either, however close the model happens to be.
    std::string dir(TEST_SAMPLES_DIR);
    const FacetedGeometry geom(dir + "/cylinder.msh");
    Blocking<FacetedGeometry> blocking(geom, 1);
    const auto b = blocking.create_hex_block(box(0, 1));

    REQUIRE(blocking.pillow(faces_of(blocking, b), b, 0.25, 1.0, 1.0, 1.0));
    REQUIRE(blocking.nb_cells<3>() == 7);
    auto &map = blocking.cmap();
    for (auto it = map.attributes<0>().begin(), end = map.attributes<0>().end(); it != end; ++it) {
        REQUIRE(it->info().geom_targets.empty());
    }
}

TEST_CASE("collapsing_a_chord_folds_its_column_away_and_joins_the_2_blocks_across_the_fold", "[BlockTestSuite]") {
    // The chord along x at (y,z) = (0,0) of a 2x2x2 grid: 2 blocks strung together through opposite
    // faces. Folding it on the diagonal through (0,0,0) takes both out, and the 2 blocks that were
    // only edge-neighbours across the fold — the one above it and the one beside it — come to share
    // a face. The valence around the chord goes from 4 to 3, which is the whole point of the
    // operation.
    const FacetedGeometry geom = make_far_away_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 1);
    build_grid_2x2x2(blocking);
    REQUIRE(blocking.nb_cells<3>() == 8);

    const auto start = face_at(blocking, {{{0, 0, 0}}, {{0, 1, 0}}, {{0, 1, 1}}, {{0, 0, 1}}});
    REQUIRE(start != nullptr);
    const auto hinge = corner_at(blocking, start, {0, 0, 0});
    REQUIRE(hinge != nullptr);

    REQUIRE(blocking.collapse_chord(start, hinge, 1e-9));
    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.nb_cells<3>() == 6);
    for (const double v : blocking.block_volumes(2)) {
        REQUIRE(v > 0.0);
    }

    // The 2 that closed the gap: the block above the chord and the one beside it.
    Blocking<FacetedGeometry>::Block above = nullptr;
    Blocking<FacetedGeometry>::Block beside = nullptr;
    auto &map = blocking.cmap();
    for (auto it = map.attributes<3>().begin(), end = map.attributes<3>().end(); it != end; ++it) {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        for (auto c = map.one_dart_per_incident_cell<0, 3>(it->dart()).begin(),
                  cend = map.one_dart_per_incident_cell<0, 3>(it->dart()).end();
             c != cend;
             ++c) {
            const Point3d &p = map.attribute<0>(c)->info().point;
            x += p.x() / 8.0;
            y += p.y() / 8.0;
            z += p.z() / 8.0;
        }
        if (x < 1.0 && y < 1.0 && z > 1.0) above = it;
        if (x < 1.0 && y > 1.0 && z < 1.0) beside = it;
    }
    REQUIRE(above != nullptr);
    REQUIRE(beside != nullptr);
    REQUIRE(share_a_face(blocking, above, beside));
}

TEST_CASE("collapsing_the_chord_of_a_lone_block_leaves_nothing", "[BlockTestSuite]") {
    // A chord of one. There is nothing either side to close over the gap, so what is left is
    // nothing — the same state `delete_sheet()` leaves when the sheet is the whole blocking.
    const FacetedGeometry geom = make_far_away_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 1);
    const auto b = blocking.create_hex_block(box(0, 1));
    const auto start = faces_of(blocking, b)[0];
    const auto hinge = blocking.cmap().attribute<0>(start->dart());

    REQUIRE(blocking.collapse_chord(start, hinge, 1e-9));
    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.nb_cells<3>() == 0);
    REQUIRE(blocking.nb_cells<0>() == 0);
}

TEST_CASE("collapsing_a_chord_refuses_what_has_no_single_fold", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_far_away_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 1);
    build_grid_2x2x2(blocking);
    const auto blocks_before = blocking.nb_cells<3>();
    const auto nodes_before = blocking.nb_cells<0>();

    const auto start = face_at(blocking, {{{0, 0, 0}}, {{0, 1, 0}}, {{0, 1, 1}}, {{0, 0, 1}}});
    REQUIRE(start != nullptr);

    SECTION("a hinge that is not a corner of the face") {
        const auto elsewhere = face_at(blocking, {{{2, 0, 0}}, {{2, 1, 0}}, {{2, 1, 1}}, {{2, 0, 1}}});
        REQUIRE(elsewhere != nullptr);
        const auto stranger = blocking.cmap().attribute<0>(elsewhere->dart());
        REQUIRE_FALSE(blocking.collapse_chord(start, stranger, 1e-9));
    }
    SECTION("a standalone quad block, which is no chord's cross-section") {
        const auto quad =
            blocking.create_quad_block({Point3d(0, 0, 5), Point3d(1, 0, 5), Point3d(1, 1, 5), Point3d(0, 1, 5)});
        const auto corner = blocking.cmap().attribute<0>(quad->dart());
        REQUIRE_FALSE(blocking.collapse_chord(quad, corner, 1e-9));
        REQUIRE(blocking.nb_cells<3>() == blocks_before);
    }
    SECTION("nothing else moved") {
        REQUIRE(blocking.nb_cells<3>() == blocks_before);
        REQUIRE(blocking.nb_cells<0>() == nodes_before);
    }
}

TEST_CASE("collapsing_a_chord_refuses_to_merge_2_different_model_vertices", "[BlockTestSuite]") {
    // The same information loss `delete_sheet()` refuses, and for the same reason: the 2 corners
    // that would meet sit on 2 different vertices of the model, and folding would leave one of them
    // with no corner of the block structure on it. Every corner of this block is on one, so every
    // fold of it is refused, whichever face and whichever diagonal it is named on.
    const FacetedGeometry geom = make_cube_vertices_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 1);
    const auto b = blocking.create_hex_block(box(0, 1));
    blocking.classify(1e-6, 1e-3, 1e-2);

    auto &map = blocking.cmap();
    int on_a_vertex = 0;
    for (auto it = map.attributes<0>().begin(), end = map.attributes<0>().end(); it != end; ++it) {
        if (!it->info().geom_targets.empty() && it->info().geom_targets[0].first == GroupDim::Dim0) {
            ++on_a_vertex;
        }
    }
    REQUIRE(on_a_vertex == 8);

    for (const auto face : faces_of(blocking, b)) {
        auto walk = face->dart();
        for (int c = 0; c < 4; ++c) {
            REQUIRE_FALSE(blocking.collapse_chord(face, map.attribute<0>(walk), 1e-6, 1e-3, 1e-2));
            walk = map.beta<1>(walk);
        }
    }
    REQUIRE(blocking.nb_cells<3>() == 1);
    REQUIRE(blocking.nb_cells<0>() == 8);
}
