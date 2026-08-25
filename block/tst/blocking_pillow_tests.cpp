#include <algorithm>
#include <array>
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
