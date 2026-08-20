#include <array>
#include <filesystem>

#include <gecko/block/Blocking.h>
#include <gecko/block/FaceCollapse.h>
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

        const auto path = (std::filesystem::temp_directory_path() / "gecko_face_collapse_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }

    /** @brief Walks AFace's dart to find the 4 corner Node handles, in local order 0..3. */
    template<typename TBlocking>
    std::array<typename TBlocking::Node, 4> corners_of(TBlocking &ABlocking, typename TBlocking::Face AFace) {
        std::array<typename TBlocking::Node, 4> corners{};
        auto d = AFace->dart();
        for (int c = 0; c < 4; ++c) {
            corners[static_cast<std::size_t>(c)] = ABlocking.cmap().template attribute<0>(d);
            d = ABlocking.cmap().template beta<1>(d);
        }
        return corners;
    }
} // namespace

TEST_CASE("collapsing_a_standalone_quad_removes_everything", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    auto face = blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
    const auto corners = corners_of(blocking, face);

    FaceCollapse<FacetedGeometry> op(blocking);
    REQUIRE(op.can_apply(face, corners[0], corners[2]));
    op.apply(face, corners[0], corners[2]);

    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.nb_cells<0>() == 0);
    REQUIRE(blocking.nb_cells<1>() == 0);
    REQUIRE(blocking.nb_cells<2>() == 0);
}

TEST_CASE("collapsing_one_of_two_sewn_quads_preserves_the_other", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    // Quad A: [0,1]x[0,1]. Quad B: [1,2]x[0,1], sharing A's edge (1,0)-(1,1) with B.
    const auto face_a = blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
    const auto face_b = blocking.create_quad_block(
        {Point3d(1.0, 0.0, 0.0), Point3d(2.0, 0.0, 0.0), Point3d(2.0, 1.0, 0.0), Point3d(1.0, 1.0, 0.0)});
    blocking.build_connectivity();
    REQUIRE(blocking.nb_cells<0>() == 6);
    REQUIRE(blocking.nb_cells<1>() == 7);
    REQUIRE(blocking.nb_cells<2>() == 2);

    const auto corners_a = corners_of(blocking, face_a);
    // corners_a walks from face_a->dart(), which starts at (0,1): corners_a =
    // {(0,1),(0,0),(1,0),(1,1)}. Collapsing the (0,1)-(1,0) diagonal keeps the shared (1,0)-(1,1)
    // edge alive through B (the other 3 boundary edges are A-only, garbage-collected).
    FaceCollapse<FacetedGeometry> op(blocking);
    REQUIRE(op.can_apply(face_a, corners_a[0], corners_a[2]));
    op.apply(face_a, corners_a[0], corners_a[2]);

    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.nb_cells<2>() == 1); // only B is left as a face
    // B itself is untouched as a quad: still exactly 4 corners, still a valid face attribute.
    const auto remaining_face = blocking.cmap().attributes<2>().begin();
    REQUIRE(remaining_face != blocking.cmap().attributes<2>().end());
    REQUIRE(remaining_face == face_b);
    int face_corner_count = 0;
    auto wd = face_b->dart();
    for (int c = 0; c < 4; ++c) {
        ++face_corner_count;
        wd = blocking.cmap().beta<1>(wd);
    }
    REQUIRE(wd == face_b->dart());
    REQUIRE(face_corner_count == 4);
}

TEST_CASE("can_apply_rejects_adjacent_corners", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    auto face = blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
    const auto corners = corners_of(blocking, face);

    FaceCollapse<FacetedGeometry> op(blocking);
    REQUIRE_FALSE(op.can_apply(face, corners[0], corners[1])); // adjacent, not diagonal
}

TEST_CASE("can_apply_rejects_when_a_hex_block_exists", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    auto face = blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
    const auto corners = corners_of(blocking, face);
    blocking.create_hex_block({Point3d(10.0, 0.0, 0.0),
                               Point3d(11.0, 0.0, 0.0),
                               Point3d(11.0, 1.0, 0.0),
                               Point3d(10.0, 1.0, 0.0),
                               Point3d(10.0, 0.0, 1.0),
                               Point3d(11.0, 0.0, 1.0),
                               Point3d(11.0, 1.0, 1.0),
                               Point3d(10.0, 1.0, 1.0)});

    FaceCollapse<FacetedGeometry> op(blocking);
    REQUIRE_FALSE(op.can_apply(face, corners[0], corners[2]));
}

TEST_CASE("can_apply_rejects_incoherent_corner_classification", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    auto face = blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
    const auto corners = corners_of(blocking, face);

    // Diagonal corners classified onto 2 genuinely different model vertices -- merging them would
    // be geometrically meaningless, so the whole collapse must be blocked.
    corners[0]->info().geom_targets = {{GroupDim::Dim0, 1}};
    corners[2]->info().geom_targets = {{GroupDim::Dim0, 2}};

    FaceCollapse<FacetedGeometry> op(blocking);
    REQUIRE_FALSE(op.can_apply(face, corners[0], corners[2]));
}

TEST_CASE("can_apply_rejects_incoherent_edge_pair_classification", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    auto face = blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
    const auto corners = corners_of(blocking, face);

    // Corners themselves are coherent (both unclassified), but the 2 edges that would merge on
    // one side of the diagonal are classified onto different curves -- still blocked.
    Blocking<FacetedGeometry>::Edge e01 = blocking.cmap().attribute<1>(face->dart());
    Blocking<FacetedGeometry>::Edge e12 = blocking.cmap().attribute<1>(blocking.cmap().beta<1>(face->dart()));
    e01->info().geom_targets = {{GroupDim::Dim1, 10}};
    e12->info().geom_targets = {{GroupDim::Dim1, 11}};

    FaceCollapse<FacetedGeometry> op(blocking);
    REQUIRE_FALSE(op.can_apply(face, corners[0], corners[2]));
}

TEST_CASE("apply_merges_node_classification_when_one_side_is_unclassified", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    // 2 sewn quads (as in collapsing_one_of_two_sewn_quads_preserves_the_other): collapsing A's
    // (0,0)-(1,1) diagonal leaves the merged node alive (via the shared edge surviving through B),
    // so there's something left to check classification on -- a fully standalone quad would
    // garbage-collect the merged node along with everything else, defeating the point of this test.
    const auto face_a = blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
    blocking.create_quad_block(
        {Point3d(1.0, 0.0, 0.0), Point3d(2.0, 0.0, 0.0), Point3d(2.0, 1.0, 0.0), Point3d(1.0, 1.0, 0.0)});
    blocking.build_connectivity();
    const auto corners = corners_of(blocking, face_a);

    corners[0]->info().geom_targets = {{GroupDim::Dim0, 42}};
    // corners[2] left unclassified (empty) -- coherent, adopts corners[0]'s classification.

    FaceCollapse<FacetedGeometry> op(blocking);
    REQUIRE(op.can_apply(face_a, corners[0], corners[2]));
    op.apply(face_a, corners[0], corners[2]);

    REQUIRE(blocking.is_valid_topology());
    bool found = false;
    for (auto it = blocking.cmap().attributes<0>().begin(), itend = blocking.cmap().attributes<0>().end(); it != itend;
         ++it) {
        if (it->info().geom_targets.size() == 1 && it->info().geom_targets[0] == std::pair{GroupDim::Dim0, Int(42)}) {
            found = true;
        }
    }
    REQUIRE(found);
}

TEST_CASE("apply_merges_edge_classification_when_one_side_is_unclassified", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    // Same 2-sewn-quads fixture. corners_of() walks from face_a->dart(), which (per this
    // CGAL version's own choice of representative dart -- not necessarily the creation-order
    // corner 0) starts at (0,1): corners = {(0,1),(0,0),(1,0),(1,1)}. Collapsing the (0,1)-(1,0)
    // diagonal makes bigon 2 = {e23=(1,0)-(1,1), e30=(1,1)-(0,1)} touch corners[3]=(1,1); e23 is
    // A's right edge, shared with B -- that's the pair with a survivor, exercised below (bigon 1
    // = {e01,e12}, both A-only, has no outside neighbor on either side and fully vanishes).
    const auto face_a = blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
    blocking.create_quad_block(
        {Point3d(1.0, 0.0, 0.0), Point3d(2.0, 0.0, 0.0), Point3d(2.0, 1.0, 0.0), Point3d(1.0, 1.0, 0.0)});
    blocking.build_connectivity();
    const auto corners = corners_of(blocking, face_a);

    auto &cmap = blocking.cmap();
    const auto d0 = face_a->dart();
    const auto d3 = cmap.beta<0>(d0);
    Blocking<FacetedGeometry>::Edge e30 = cmap.attribute<1>(d3);
    e30->info().geom_targets = {{GroupDim::Dim1, 99}};
    // e23 (the other bigon-2 edge, shared with B) left unclassified -- coherent, adopts e30's tag.

    FaceCollapse<FacetedGeometry> op(blocking);
    REQUIRE(op.can_apply(face_a, corners[0], corners[2]));
    op.apply(face_a, corners[0], corners[2]);

    REQUIRE(blocking.is_valid_topology());
    bool found = false;
    for (auto it = blocking.cmap().attributes<1>().begin(), itend = blocking.cmap().attributes<1>().end(); it != itend;
         ++it) {
        if (it->info().geom_targets.size() == 1 && it->info().geom_targets[0] == std::pair{GroupDim::Dim1, Int(99)}) {
            found = true;
        }
    }
    REQUIRE(found);
}
