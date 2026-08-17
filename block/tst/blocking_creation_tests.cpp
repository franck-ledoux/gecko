#include <array>
#include <filesystem>
#include <set>
#include <utility>

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

        const auto path = (std::filesystem::temp_directory_path() / "gecko_block_creation_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }

    /**
     * @brief Collects every edge of @p ABlocking as a (min,max) pair of indices into @p ACorners
     * (matched by exact position), for order-independent comparison against an expected edge set.
     */
    template<typename TBlocking>
    std::set<std::pair<int, int>> edge_index_pairs(TBlocking &ABlocking, const std::array<Point3d, 8> &ACorners) {
        auto index_of = [&](const Point3d &p) -> int {
            for (int i = 0; i < 8; ++i) {
                if (ACorners[static_cast<std::size_t>(i)] == p) return i;
            }
            FAIL("corner position not found among ACorners");
            return -1;
        };

        std::set<std::pair<int, int>> pairs;
        for (auto it = ABlocking.cmap().template attributes<1>().begin(),
                  itend = ABlocking.cmap().template attributes<1>().end();
             it != itend;
             ++it) {
            auto d = it->dart();
            const int a = index_of(ABlocking.cmap().template attribute<0>(d)->info().point);
            const int b =
                index_of(ABlocking.cmap().template attribute<0>(ABlocking.cmap().template beta<1>(d))->info().point);
            pairs.insert(a < b ? std::make_pair(a, b) : std::make_pair(b, a));
        }
        return pairs;
    }

    /**
     * @brief Collects every face of @p ABlocking as a sorted set of indices into @p ACorners
     * (matched by exact position), for order-independent comparison against expected face sets.
     */
    template<typename TBlocking>
    std::set<std::set<int>> face_index_sets(TBlocking &ABlocking, const std::array<Point3d, 8> &ACorners) {
        auto index_of = [&](const Point3d &p) -> int {
            for (int i = 0; i < 8; ++i) {
                if (ACorners[static_cast<std::size_t>(i)] == p) return i;
            }
            FAIL("corner position not found among ACorners");
            return -1;
        };

        std::set<std::set<int>> faces;
        for (auto it = ABlocking.cmap().template attributes<2>().begin(),
                  itend = ABlocking.cmap().template attributes<2>().end();
             it != itend;
             ++it) {
            auto d = it->dart();
            std::set<int> face;
            for (int k = 0; k < 4; ++k) {
                face.insert(index_of(ABlocking.cmap().template attribute<0>(d)->info().point));
                d = ABlocking.cmap().template beta<1>(d);
            }
            faces.insert(face);
        }
        return faces;
    }
} // namespace

TEST_CASE("quad_block_creation_has_no_3_cell", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    const std::array<Point3d, 4> corners = {
        Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)};
    auto face = blocking.create_quad_block(corners);

    REQUIRE(blocking.nb_cells<0>() == 4);
    REQUIRE(blocking.nb_cells<1>() == 4);
    REQUIRE(blocking.nb_cells<2>() == 1);
    REQUIRE(blocking.nb_cells<3>() == 0); // a "2D block" is a standalone face, no 3-cell created
    REQUIRE(blocking.is_valid_topology());

    // The face's own Coons surface (N=1, so bilinear) must reproduce all 4 corners and the center.
    REQUIRE(face->info().surface.value(0.0, 0.0) == corners[0]);
    REQUIRE(face->info().surface.value(1.0, 0.0) == corners[1]);
    REQUIRE(face->info().surface.value(1.0, 1.0) == corners[2]);
    REQUIRE(face->info().surface.value(0.0, 1.0) == corners[3]);
    const Point3d center = face->info().surface.value(0.5, 0.5);
    REQUIRE(center.x() == Approx(0.5));
    REQUIRE(center.y() == Approx(0.5));
    REQUIRE(center.z() == Approx(0.0));
}

TEST_CASE("hex_block_creation_has_expected_cell_counts", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    const std::array<Point3d, 8> corners = {
        Point3d(0.0, 0.0, 0.0),
        Point3d(1.0, 0.0, 0.0),
        Point3d(1.0, 1.0, 0.0),
        Point3d(0.0, 1.0, 0.0), // bottom
        Point3d(0.0, 0.0, 1.0),
        Point3d(1.0, 0.0, 1.0),
        Point3d(1.0, 1.0, 1.0),
        Point3d(0.0, 1.0, 1.0) // top
    };
    blocking.create_hex_block(corners);

    REQUIRE(blocking.nb_cells<0>() == 8);
    REQUIRE(blocking.nb_cells<1>() == 12);
    REQUIRE(blocking.nb_cells<2>() == 6);
    REQUIRE(blocking.nb_cells<3>() == 1);
    REQUIRE(blocking.is_valid_topology());
}

TEST_CASE("hex_block_corner_order_matches_gecko_hex8_edge_convention", "[BlockTestSuite]") {
    // Empirically verifies (rather than assumes) that create_hex_block()'s corner-index-to-node
    // correspondence, for corners given in gecko::CubicTraits::Cell's own documented HEX8 order
    // (bottom perimeter 0-1-2-3, top perimeter 4-5-6-7 directly above), actually produces the same
    // 12 edges CubicTraits::Cell::EdgeNodes documents — i.e. block/'s internal convention happens to
    // coincide with the rest of gecko's, which to_mesh() (a later task) can then rely on directly.
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    const std::array<Point3d, 8> corners = {
        Point3d(0.0, 0.0, 0.0),
        Point3d(1.0, 0.0, 0.0),
        Point3d(1.0, 1.0, 0.0),
        Point3d(0.0, 1.0, 0.0), // bottom
        Point3d(0.0, 0.0, 1.0),
        Point3d(1.0, 0.0, 1.0),
        Point3d(1.0, 1.0, 1.0),
        Point3d(0.0, 1.0, 1.0) // top
    };
    blocking.create_hex_block(corners);

    const std::set<std::pair<int, int>> actual_edges = edge_index_pairs(blocking, corners);
    const std::set<std::pair<int, int>> expected_edges = {
        {0, 1},
        {1, 2},
        {2, 3},
        {0, 3}, // bottom perimeter
        {4, 5},
        {5, 6},
        {6, 7},
        {4, 7}, // top perimeter
        {0, 4},
        {1, 5},
        {2, 6},
        {3, 7} // vertical
    };
    REQUIRE(actual_edges == expected_edges);

    // Same check for face grouping, against CubicTraits::Cell::FaceNodes (order-independent: only
    // which 4 corners group into a face is checked here, not their cyclic order/orientation).
    const std::set<std::set<int>> actual_faces = face_index_sets(blocking, corners);
    const std::set<std::set<int>> expected_faces = {
        {0, 1, 2, 3}, // bottom
        {4, 5, 6, 7}, // top
        {0, 1, 4, 5}, // front
        {1, 2, 5, 6}, // right
        {2, 3, 6, 7}, // back
        {0, 3, 4, 7}  // left
    };
    REQUIRE(actual_faces == expected_faces);
}

TEST_CASE("hex_block_volume_and_faces_reproduce_straight_geometry", "[BlockTestSuite]") {
    // End-to-end check of the HEX_EDGES/HEX_FACES orientation tables in Blocking.h: for a unit-cube
    // block (N=1, so Coons/TFI collapse to bilinear/trilinear interpolation), the TFI-evaluated
    // volume must reproduce known corners and the geometric center exactly. A wrong table entry
    // (bad node index, wrong Fu/Fv/Fw role, or a direction that needed reversing) would make this
    // fail, even though the cell-count/topology checks above would still pass.
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    const std::array<Point3d, 8> corners = {
        Point3d(0.0, 0.0, 0.0),
        Point3d(1.0, 0.0, 0.0),
        Point3d(1.0, 1.0, 0.0),
        Point3d(0.0, 1.0, 0.0), // bottom
        Point3d(0.0, 0.0, 1.0),
        Point3d(1.0, 0.0, 1.0),
        Point3d(1.0, 1.0, 1.0),
        Point3d(0.0, 1.0, 1.0) // top
    };
    auto block = blocking.create_hex_block(corners);

    // (u,v,w) = (0,0,0) is the corner common to Fu0/Fv0/Fw0 -> node 0; (1,1,1) common to
    // Fu1/Fv1/Fw1 -> node 6; (1,0,0) -> node 1; (0,1,0) -> node 3 (derived by hand from HEX_FACES).
    REQUIRE(block->info().volume.value(0.0, 0.0, 0.0) == corners[0]);
    REQUIRE(block->info().volume.value(1.0, 1.0, 1.0) == corners[6]);
    REQUIRE(block->info().volume.value(1.0, 0.0, 0.0) == corners[1]);
    REQUIRE(block->info().volume.value(0.0, 1.0, 0.0) == corners[3]);

    const Point3d center = block->info().volume.value(0.5, 0.5, 0.5);
    REQUIRE(center.x() == Approx(0.5));
    REQUIRE(center.y() == Approx(0.5));
    REQUIRE(center.z() == Approx(0.5));

    // Every one of the 6 face surfaces must also independently reproduce its own corners exactly.
    for (auto it = blocking.cmap().template attributes<2>().begin(),
              itend = blocking.cmap().template attributes<2>().end();
         it != itend;
         ++it) {
        const auto &surface = it->info().surface;
        const Point3d p00 = surface.value(0.0, 0.0);
        const Point3d p11 = surface.value(1.0, 1.0);
        // Every face corner value must be one of the 8 known cube corners (exact match).
        bool p00_is_known = false, p11_is_known = false;
        for (const auto &c : corners) {
            if (p00 == c) p00_is_known = true;
            if (p11 == c) p11_is_known = true;
        }
        REQUIRE(p00_is_known);
        REQUIRE(p11_is_known);
    }
}
