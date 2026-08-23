#include <array>
#include <cmath>
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

        const auto path = (std::filesystem::temp_directory_path() / "gecko_block_connectivity_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }
} // namespace

TEST_CASE("two_hex_blocks_sew_at_shared_face", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    // Block A: unit cube [0,1]^3. Block B: unit cube [1,2]x[0,1]x[0,1], sharing A's x=1 face
    // (A's Fu1 = corners 1,2,5,6) with B's x=1 face (B's Fu0 = corners 0,3,4,7).
    const std::array<Point3d, 8> corners_a = {Point3d(0.0, 0.0, 0.0),
                                              Point3d(1.0, 0.0, 0.0),
                                              Point3d(1.0, 1.0, 0.0),
                                              Point3d(0.0, 1.0, 0.0),
                                              Point3d(0.0, 0.0, 1.0),
                                              Point3d(1.0, 0.0, 1.0),
                                              Point3d(1.0, 1.0, 1.0),
                                              Point3d(0.0, 1.0, 1.0)};
    const std::array<Point3d, 8> corners_b = {Point3d(1.0, 0.0, 0.0),
                                              Point3d(2.0, 0.0, 0.0),
                                              Point3d(2.0, 1.0, 0.0),
                                              Point3d(1.0, 1.0, 0.0),
                                              Point3d(1.0, 0.0, 1.0),
                                              Point3d(2.0, 0.0, 1.0),
                                              Point3d(2.0, 1.0, 1.0),
                                              Point3d(1.0, 1.0, 1.0)};

    blocking.create_hex_block(corners_a);
    blocking.create_hex_block(corners_b);

    // Before sewing: 16 nodes (8+8, none shared yet), 2 separate blocks.
    REQUIRE(blocking.nb_cells<0>() == 16);
    REQUIRE(blocking.nb_cells<3>() == 2);

    blocking.build_connectivity();

    REQUIRE(blocking.is_valid_topology());
    // After sewing at the shared face: the 4 shared corners are merged (16 - 4 = 12 nodes), the 4
    // shared edges are merged (12+12-4=20 edges), and the 1 shared face is merged (6+6-1=11 faces);
    // the 2 blocks themselves remain distinct 3-cells.
    REQUIRE(blocking.nb_cells<0>() == 12);
    REQUIRE(blocking.nb_cells<1>() == 20);
    REQUIRE(blocking.nb_cells<2>() == 11);
    REQUIRE(blocking.nb_cells<3>() == 2);

    // Stronger than a count match: a shared corner (e.g. (1,0,0), present in both blocks' original
    // corner arrays) must now be exactly ONE node attribute, not two coincidentally-equal ones.
    int count_at_shared_corner = 0;
    for (auto it = blocking.cmap().attributes<0>().begin(), itend = blocking.cmap().attributes<0>().end(); it != itend;
         ++it) {
        if (it->info().point == Point3d(1.0, 0.0, 0.0)) ++count_at_shared_corner;
    }
    REQUIRE(count_at_shared_corner == 1);
}

TEST_CASE("two_quad_blocks_sew_at_shared_edge", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    // Quad A: [0,1]x[0,1]. Quad B: [1,2]x[0,1], sharing A's edge (corner1->corner2) with B's edge
    // (corner3->corner0).
    const std::array<Point3d, 4> corners_a = {
        Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)};
    const std::array<Point3d, 4> corners_b = {
        Point3d(1.0, 0.0, 0.0), Point3d(2.0, 0.0, 0.0), Point3d(2.0, 1.0, 0.0), Point3d(1.0, 1.0, 0.0)};

    blocking.create_quad_block(corners_a);
    blocking.create_quad_block(corners_b);

    REQUIRE(blocking.nb_cells<0>() == 8);
    REQUIRE(blocking.nb_cells<2>() == 2);

    blocking.build_connectivity();

    REQUIRE(blocking.is_valid_topology());
    // 2 shared corners merged (8-2=6 nodes), 1 shared edge merged (4+4-1=7 edges), 2 faces remain,
    // still no 3-cells (a "2D block structure" stays 2D).
    REQUIRE(blocking.nb_cells<0>() == 6);
    REQUIRE(blocking.nb_cells<1>() == 7);
    REQUIRE(blocking.nb_cells<2>() == 2);
    REQUIRE(blocking.nb_cells<3>() == 0);
}

TEST_CASE("non_adjacent_blocks_stay_unsewn", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    const std::array<Point3d, 8> corners_a = {Point3d(0.0, 0.0, 0.0),
                                              Point3d(1.0, 0.0, 0.0),
                                              Point3d(1.0, 1.0, 0.0),
                                              Point3d(0.0, 1.0, 0.0),
                                              Point3d(0.0, 0.0, 1.0),
                                              Point3d(1.0, 0.0, 1.0),
                                              Point3d(1.0, 1.0, 1.0),
                                              Point3d(0.0, 1.0, 1.0)};
    // Far away, no shared geometry at all.
    const std::array<Point3d, 8> corners_b = {Point3d(10.0, 0.0, 0.0),
                                              Point3d(11.0, 0.0, 0.0),
                                              Point3d(11.0, 1.0, 0.0),
                                              Point3d(10.0, 1.0, 0.0),
                                              Point3d(10.0, 0.0, 1.0),
                                              Point3d(11.0, 0.0, 1.0),
                                              Point3d(11.0, 1.0, 1.0),
                                              Point3d(10.0, 1.0, 1.0)};

    blocking.create_hex_block(corners_a);
    blocking.create_hex_block(corners_b);
    blocking.build_connectivity();

    REQUIRE(blocking.is_valid_topology());
    // Nothing shared: counts stay exactly the sum of 2 independent hexes.
    REQUIRE(blocking.nb_cells<0>() == 16);
    REQUIRE(blocking.nb_cells<1>() == 24);
    REQUIRE(blocking.nb_cells<2>() == 12);
    REQUIRE(blocking.nb_cells<3>() == 2);
}

TEST_CASE("sewing_glues_2_blocks_corner_to_corner_rather_than_a_quarter_turn_out", "[BlockTestSuite]") {
    // Sewing links darts, and the 2 it links are opposite halves of one edge — so the dart matching
    // a given one is the dart whose own start is that one's *end*. Testing the alignment against the
    // dart's start instead is one step out of phase, and because the search then rotates its
    // candidate until something matches, it does not fail: it succeeds on the neighbouring rotation
    // and glues every shared face a quarter turn twisted.
    //
    // Nothing noticed for a long time. A block's volume is built once from the corners it was created
    // with, and to_mesh() reconciles the 2 sides of a face by position, so almost nothing reads a
    // block's structure back out of its darts. Anything that does — move_node(), delete_sheet() —
    // gets a block whose opposite faces disagree about which way v and w run.
    //
    // Stated here as the property that catches it whatever the cause: on an axis-aligned blocking,
    // the 2 corners of every edge differ in exactly 1 coordinate.
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_hex_block({Point3d(0, 0, 0),
                               Point3d(1, 0, 0),
                               Point3d(1, 1, 0),
                               Point3d(0, 1, 0),
                               Point3d(0, 0, 1),
                               Point3d(1, 0, 1),
                               Point3d(1, 1, 1),
                               Point3d(0, 1, 1)});
    blocking.create_hex_block({Point3d(1, 0, 0),
                               Point3d(2, 0, 0),
                               Point3d(2, 1, 0),
                               Point3d(1, 1, 0),
                               Point3d(1, 0, 1),
                               Point3d(2, 0, 1),
                               Point3d(2, 1, 1),
                               Point3d(1, 1, 1)});

    auto &map = blocking.cmap();
    const auto skew_darts = [&map]() {
        int skew = 0;
        for (auto d = map.darts().begin(), end = map.darts().end(); d != end; ++d) {
            const Point3d a = map.attribute<0>(d)->info().point;
            const Point3d b = map.attribute<0>(map.beta<1>(d))->info().point;
            int differing = 0;
            if (std::abs(a.x() - b.x()) > 1e-9) ++differing;
            if (std::abs(a.y() - b.y()) > 1e-9) ++differing;
            if (std::abs(a.z() - b.z()) > 1e-9) ++differing;
            if (differing != 1) ++skew;
        }
        return skew;
    };

    REQUIRE(skew_darts() == 0);
    blocking.build_connectivity();
    REQUIRE(skew_darts() == 0);
    REQUIRE(blocking.is_valid_topology());

    // And the consequence the property is there to protect: a block rebuilt from its own darts still
    // measures what it looks like. move_node() to where the corner already is rebuilds every cell
    // touching it without moving anything.
    const auto corner = map.attributes<0>().begin();
    blocking.move_node(corner, corner->info().point);
    for (const double v : blocking.block_volumes(2)) {
        REQUIRE(v > 0.99);
        REQUIRE(v < 1.01);
    }
}
