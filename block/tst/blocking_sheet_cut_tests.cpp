#include <array>
#include <cmath>
#include <filesystem>
#include <vector>

#include <gecko/block/Blocking.h>
#include <gecko/geom/FacetedGeometry.h>
#include <gecko/io/GmshMeshWriter.h>
#include <gecko/math/Vector3d.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
using Catch::Approx;

using namespace gecko;

namespace {
    /** @brief Builds a minimal FacetedGeometry fixture: a single tagged triangle, deliberately far
     * from every block below so nothing ever classifies and the cut is exercised on its own. */
    FacetedGeometry make_far_geom_model() {
        SimplicialMesh mesh;
        auto n0 = mesh.add_node(100.0, 100.0, 100.0);
        auto n1 = mesh.add_node(101.0, 100.0, 100.0);
        auto n2 = mesh.add_node(100.0, 101.0, 100.0);

        GroupRegistry groups;
        auto surf = groups.add_group("Surf", GroupDim::Dim2);
        auto &face_group = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &face_entity = mesh.add_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
        auto f0 = mesh.add_face(n0, n1, n2);
        face_group[f0.value] = surf;
        face_entity[f0.value] = 1;

        const auto path = (std::filesystem::temp_directory_path() / "gecko_sheet_cut_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }

    /** @brief The 8 corners of an axis-aligned box, in `create_hex_block`'s expected order. */
    std::array<Point3d, 8> box(double AX0, double AX1) {
        return {Point3d(AX0, 0.0, 0.0),
                Point3d(AX1, 0.0, 0.0),
                Point3d(AX1, 1.0, 0.0),
                Point3d(AX0, 1.0, 0.0),
                Point3d(AX0, 0.0, 1.0),
                Point3d(AX1, 0.0, 1.0),
                Point3d(AX1, 1.0, 1.0),
                Point3d(AX0, 1.0, 1.0)};
    }

    /** @brief The first edge of @p ABlocking whose 2 endpoints differ only along @p AAxis — i.e. an
     * edge running along that world axis, used to aim a cut without depending on traversal order. */
    template<typename TBlocking>
    typename TBlocking::Edge edge_along(TBlocking &ABlocking, int AAxis, double AOtherCoordSum) {
        auto &map = ABlocking.cmap();
        for (auto it = map.template attributes<1>().begin(), itend = map.template attributes<1>().end(); it != itend;
             ++it) {
            const auto d = it->dart();
            const Point3d &a = map.template attribute<0>(d)->info().point;
            const Point3d &b = map.template attribute<0>(map.template beta<1>(d))->info().point;
            const std::array<double, 3> pa{a.x(), a.y(), a.z()};
            const std::array<double, 3> pb{b.x(), b.y(), b.z()};
            bool runs_along = true;
            double other_sum = 0.0;
            for (int k = 0; k < 3; ++k) {
                const bool same = std::abs(pa[k] - pb[k]) < 1e-12;
                if (k == AAxis && same) runs_along = false;
                if (k != AAxis) {
                    if (!same) runs_along = false;
                    other_sum += pa[k];
                }
            }
            if (runs_along && std::abs(other_sum - AOtherCoordSum) < 1e-12) return it;
        }
        FAIL("no edge running along the requested axis");
        return {};
    }
} // namespace

TEST_CASE("sheet_through_a_lone_hex_edge_is_its_4_parallel_edges", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_hex_block(box(0.0, 1.0));

    // An x-running edge: its sheet is the block's 4 x-edges, the 4 faces holding them (every face
    // except x=0 and x=1), and the block itself.
    const auto e = edge_along(blocking, 0, 0.0);
    const auto sheet = blocking.find_sheet(e);

    REQUIRE(sheet.has_value());
    REQUIRE(sheet->edges.size() == 4);
    REQUIRE(sheet->faces.size() == 4);
    REQUIRE(sheet->blocks.size() == 1);
}

TEST_CASE("sheet_propagates_into_the_neighbour_block_only_across_the_shared_face", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_hex_block(box(0.0, 1.0));
    blocking.create_hex_block(box(1.0, 2.0));
    blocking.build_connectivity();
    REQUIRE(blocking.is_valid_topology());

    SECTION("a y-running edge drags both blocks in, through the 2 y-edges they share") {
        // The 2 blocks share the x=1 face, whose 4 edges include 2 y-running ones. Those are single
        // Edge attributes belonging to both blocks, so the sheet crosses over: 4 + 4 - 2 = 6 edges,
        // 4 + 4 - 1 = 7 faces (the shared one counted once), 2 blocks.
        const auto e = edge_along(blocking, 1, 0.0);
        const auto sheet = blocking.find_sheet(e);

        REQUIRE(sheet.has_value());
        REQUIRE(sheet->edges.size() == 6);
        REQUIRE(sheet->faces.size() == 7);
        REQUIRE(sheet->blocks.size() == 2);
    }

    SECTION("an x-running edge stops at the shared face, which it never crosses") {
        // The x-edges of the 2 blocks are disjoint (x in [0,1] vs [1,2]), so nothing propagates —
        // cutting one block along x leaves the shared x=1 face whole, and stays conformal.
        const auto e = edge_along(blocking, 0, 0.0);
        const auto sheet = blocking.find_sheet(e);

        REQUIRE(sheet.has_value());
        REQUIRE(sheet->edges.size() == 4);
        REQUIRE(sheet->faces.size() == 4);
        REQUIRE(sheet->blocks.size() == 1);
    }
}

TEST_CASE("sheet_through_a_quad_blocking_is_the_2D_chord", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_quad_block(
        {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)});
    blocking.create_quad_block(
        {Point3d(1.0, 0.0, 0.0), Point3d(2.0, 0.0, 0.0), Point3d(2.0, 1.0, 0.0), Point3d(1.0, 1.0, 0.0)});
    blocking.build_connectivity();
    REQUIRE(blocking.is_valid_topology());

    SECTION("a y-running edge crosses the shared edge into both quads") {
        // The 2 quads share their x=1 edge, which runs along y: 2 + 2 - 1 = 3 edges, both faces.
        const auto e = edge_along(blocking, 1, 0.0);
        const auto sheet = blocking.find_sheet(e);

        REQUIRE(sheet.has_value());
        REQUIRE(sheet->edges.size() == 3);
        REQUIRE(sheet->faces.size() == 2);
        REQUIRE(sheet->blocks.empty());
    }

    SECTION("an x-running edge stays in its own quad") {
        const auto e = edge_along(blocking, 0, 0.0);
        const auto sheet = blocking.find_sheet(e);

        REQUIRE(sheet.has_value());
        REQUIRE(sheet->edges.size() == 2);
        REQUIRE(sheet->faces.size() == 1);
        REQUIRE(sheet->blocks.empty());
    }
}

TEST_CASE("cutting_a_lone_hex_gives_2_conformal_hexes", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_hex_block(box(0.0, 1.0));

    const auto e = edge_along(blocking, 0, 0.0);
    REQUIRE(blocking.cut_sheet(e, 0.25));

    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.template nb_cells<3>() == 2);
    REQUIRE(blocking.template nb_cells<0>() == 12); // 8 corners + 4 cut nodes
    REQUIRE(blocking.template nb_cells<1>() == 20); // 12 edges, 4 of them split, + 4 new
    REQUIRE(blocking.template nb_cells<2>() == 11); // 6 faces, 4 of them split, + 1 new
}

TEST_CASE("cut_sheet_rejects_out_of_range_parameters_without_touching_anything", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_hex_block(box(0.0, 1.0));
    const auto e = edge_along(blocking, 0, 0.0);

    REQUIRE_FALSE(blocking.cut_sheet(e, 0.0));
    REQUIRE_FALSE(blocking.cut_sheet(e, 1.0));
    REQUIRE_FALSE(blocking.cut_sheet(e, -0.5));
    REQUIRE(blocking.template nb_cells<3>() == 1);
    REQUIRE(blocking.template nb_cells<0>() == 8);
}

namespace {
    /** @brief A deliberately non-affine hex: no face is planar and the trilinear map is genuinely
     * nonlinear, so a cut that picked the wrong axis, the wrong side or the wrong frame lands its
     * nodes somewhere the parent block never passes through. An axis-aligned box could not tell
     * those apart — its parameter grid is a regular spatial grid whichever way you slice it. */
    std::array<Point3d, 8> twisted_box() {
        return {Point3d(0.0, 0.0, 0.0),
                Point3d(1.0, 0.0, 0.0),
                Point3d(1.2, 1.1, 0.1),
                Point3d(0.0, 1.0, 0.0),
                Point3d(0.0, 0.0, 1.0),
                Point3d(1.0, 0.15, 1.1),
                Point3d(1.3, 1.4, 1.2),
                Point3d(0.1, 1.0, 1.0)};
    }

    /** @brief The edge of @p ABlocking joining the 2 given positions, in either order. */
    template<typename TBlocking>
    typename TBlocking::Edge edge_between(TBlocking &ABlocking, const Point3d &AA, const Point3d &AB) {
        auto &map = ABlocking.cmap();
        for (auto it = map.template attributes<1>().begin(), itend = map.template attributes<1>().end(); it != itend;
             ++it) {
            const auto d = it->dart();
            const Point3d &p = map.template attribute<0>(d)->info().point;
            const Point3d &q = map.template attribute<0>(map.template beta<1>(d))->info().point;
            if ((p == AA && q == AB) || (p == AB && q == AA)) return it;
        }
        FAIL("no edge between the 2 requested positions");
        return {};
    }

    /** @brief Whether @p APoint coincides with one of @p ACloud, within @p ATol. */
    bool lies_in(const Point3d &APoint, const std::vector<Point3d> &ACloud, double ATol) {
        for (const Point3d &c : ACloud) {
            if (Vector3d(APoint, c).norm() <= ATol) return true;
        }
        return false;
    }
} // namespace

TEST_CASE("cutting_leaves_every_generated_node_on_the_uncut_geometry", "[BlockTestSuite]") {
    // The invariant the whole operation rests on: a cut re-parameterizes the blocking without moving
    // it. Cutting at 1/3 and meshing each half with 4 intervals samples the parent at parameters
    // 0,1/12,...,4/12 on the near side and 4/12,...,12/12 on the far one, and at 0,3/12,...,12/12
    // across the 2 axes the cut left alone — every one of them a node the uncut blocking's own
    // 12-interval mesh already carries. So the cut mesh must be a strict subset of the uncut one.
    const FacetedGeometry geom = make_far_geom_model();

    Blocking<FacetedGeometry> reference(geom);
    reference.create_hex_block(twisted_box());
    const auto fine = reference.to_mesh(12);
    std::vector<Point3d> cloud;
    cloud.reserve(fine.nb_nodes());
    for (UInt i = 0; i < fine.nb_nodes(); ++i) {
        cloud.push_back(fine.node(NodeId{i}));
    }

    Blocking<FacetedGeometry> cut(geom);
    cut.create_hex_block(twisted_box());
    const auto e = edge_between(cut, Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0));
    REQUIRE(cut.cut_sheet(e, 1.0 / 3.0));
    REQUIRE(cut.is_valid_topology());

    const auto coarse = cut.to_mesh(4);
    REQUIRE(coarse.nb_cells() == 2 * 4 * 4 * 4);
    for (UInt i = 0; i < coarse.nb_nodes(); ++i) {
        REQUIRE(lies_in(coarse.node(NodeId{i}), cloud, 1e-9));
    }
}

TEST_CASE("cutting_a_curved_blocking_keeps_the_curvature_rather_than_refitting_it", "[BlockTestSuite]") {
    // The case that separates subdividing the stored geometry from rebuilding it through
    // coons_surface_from_edges()/tfi_volume_from_faces(): on a straight-edged block the 2 agree
    // exactly, so only a genuinely curved one can tell them apart. A refit would pull the interior
    // cut face onto the Coons patch of its 4 new boundary curves, which is not the surface the
    // parent block actually carried there, and its nodes would leave the parent's geometry.
    //
    // The 3 stored geometries are bent independently rather than kept mutually consistent, because
    // that is exactly how the cut treats them: an edge's curve, a face's surface and a block's
    // volume are each subdivided from their own stored form, never re-derived from one another.
    // Only strictly interior control points move, so every cell still meets its corners where the
    // Node attributes say it does.
    const FacetedGeometry geom = make_far_geom_model();
    using CubicBlocking = Blocking<FacetedGeometry, BezierCurve<3, Point3d>>;

    // Driven purely by each control point's own position, so the 2 blockings below come out
    // identical without depending on any traversal order.
    const auto wobble = [](const Point3d &AP) {
        return AP + Vector3d(0.18 * std::sin(3.0 * AP.y() + 1.0),
                             0.15 * std::sin(3.0 * AP.z() + 2.0),
                             0.20 * std::sin(3.0 * AP.x() + 0.5));
    };
    const auto bend = [&wobble](CubicBlocking &ABlocking) {
        auto &map = ABlocking.cmap();
        for (auto it = map.attributes<1>().begin(), itend = map.attributes<1>().end(); it != itend; ++it) {
            for (std::size_t i = 1; i < 3; ++i) {
                it->info().curve[i] = wobble(it->info().curve[i]);
            }
        }
        for (auto it = map.attributes<2>().begin(), itend = map.attributes<2>().end(); it != itend; ++it) {
            for (std::size_t i = 1; i < 3; ++i) {
                for (std::size_t j = 1; j < 3; ++j) {
                    it->info().surface.control_point(i, j) = wobble(it->info().surface.control_point(i, j));
                }
            }
        }
        for (auto it = map.attributes<3>().begin(), itend = map.attributes<3>().end(); it != itend; ++it) {
            for (std::size_t i = 1; i < 3; ++i) {
                for (std::size_t j = 1; j < 3; ++j) {
                    for (std::size_t k = 1; k < 3; ++k) {
                        it->info().volume.control_point(i, j, k) = wobble(it->info().volume.control_point(i, j, k));
                    }
                }
            }
        }
    };

    CubicBlocking reference(geom);
    reference.create_hex_block(twisted_box());
    bend(reference);

    CubicBlocking cut(geom);
    cut.create_hex_block(twisted_box());
    bend(cut);

    // Guard against silently testing a straight block: the bend has to reach the sampled geometry,
    // not just sit in control points the evaluation never notices.
    CubicBlocking unbent(geom);
    unbent.create_hex_block(twisted_box());
    const auto bent_mesh = reference.to_mesh(3);
    const auto straight_mesh = unbent.to_mesh(3);
    REQUIRE(bent_mesh.nb_nodes() == straight_mesh.nb_nodes());
    double largest_shift = 0.0;
    for (UInt i = 0; i < bent_mesh.nb_nodes(); ++i) {
        largest_shift =
            std::max(largest_shift, Vector3d(bent_mesh.node(NodeId{i}), straight_mesh.node(NodeId{i})).norm());
    }
    REQUIRE(largest_shift > 1e-2);

    const auto fine = reference.to_mesh(12);
    std::vector<Point3d> cloud;
    cloud.reserve(fine.nb_nodes());
    for (UInt i = 0; i < fine.nb_nodes(); ++i) {
        cloud.push_back(fine.node(NodeId{i}));
    }

    const auto e = edge_between(cut, Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0));
    REQUIRE(cut.cut_sheet(e, 1.0 / 3.0));
    REQUIRE(cut.is_valid_topology());

    const auto coarse = cut.to_mesh(4);
    for (UInt i = 0; i < coarse.nb_nodes(); ++i) {
        REQUIRE(lies_in(coarse.node(NodeId{i}), cloud, 1e-9));
    }
}
TEST_CASE("cutting_2_sewn_hexes_splits_both_and_stays_conformal", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_hex_block(box(0.0, 1.0));
    blocking.create_hex_block(box(1.0, 2.0));
    blocking.build_connectivity();

    // A y-running edge: the sheet reaches both blocks (6 edges, 7 faces, 2 blocks), so the cut has
    // to split both or the shared face would come apart.
    const auto e = edge_along(blocking, 1, 0.0);
    REQUIRE(blocking.cut_sheet(e, 0.5));

    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.template nb_cells<3>() == 4);
    // 12 corners + 1 node per cut edge.
    REQUIRE(blocking.template nb_cells<0>() == 12 + 6);
    // The 2 halves of every block still meet the neighbour's, so the mesh stays conformal: 4 blocks
    // at 1 interval share their faces rather than duplicating nodes.
    const auto mesh = blocking.to_mesh(1);
    REQUIRE(mesh.nb_cells() == 4);
    REQUIRE(mesh.nb_nodes() == 18);
}

namespace {
    /** @brief The 1x1x2 block these last tests exercise, in `create_hex_block`'s expected HEX8
     * order: extreme points (0,0,0) and (1,1,2), long axis z. */
    std::array<Point3d, 8> tall_box() {
        return {Point3d(0.0, 0.0, 0.0),
                Point3d(1.0, 0.0, 0.0),
                Point3d(1.0, 1.0, 0.0),
                Point3d(0.0, 1.0, 0.0),
                Point3d(0.0, 0.0, 2.0),
                Point3d(1.0, 0.0, 2.0),
                Point3d(1.0, 1.0, 2.0),
                Point3d(0.0, 1.0, 2.0)};
    }
} // namespace

TEST_CASE("cutting_a_1x1x2_block_in_half_gives_2_blocks_of_volume_1", "[BlockTestSuite]") {
    // Volume is what makes this bite: the corners of a half can be right while the parametrization
    // between them is wrong, and only integrating the block's own stored geometry notices.
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_hex_block(tall_box());
    REQUIRE(blocking.block_volumes()[0] == Approx(2.0).margin(1e-12));

    const auto e = edge_between(blocking, Point3d(0.0, 0.0, 0.0), Point3d(0.0, 0.0, 2.0));
    REQUIRE(blocking.cut_sheet(e, 0.5));

    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.template nb_cells<3>() == 2);

    // Sampled finely as well as coarsely: a half carrying the wrong sub-volume, or the right one in
    // a rotated frame, still agrees with its corners and disagrees everywhere between them.
    for (const SizeT s : {SizeT(1), SizeT(2), SizeT(4)}) {
        const auto volumes = blocking.block_volumes(s);
        REQUIRE(volumes.size() == 2);
        REQUIRE(volumes[0] == Approx(1.0).margin(1e-12));
        REQUIRE(volumes[1] == Approx(1.0).margin(1e-12));
    }

    // The 2 halves must meet along one face carrying exactly the 4 midplane corners.
    std::vector<Point3d> midplane;
    for (auto it = blocking.cmap().attributes<0>().begin(), end = blocking.cmap().attributes<0>().end(); it != end;
         ++it) {
        if (std::abs(it->info().point.z() - 1.0) < 1e-12) midplane.push_back(it->info().point);
    }
    REQUIRE(midplane.size() == 4);
    for (const Point3d &want :
         {Point3d(0.0, 0.0, 1.0), Point3d(0.0, 1.0, 1.0), Point3d(1.0, 1.0, 1.0), Point3d(1.0, 0.0, 1.0)}) {
        REQUIRE(lies_in(want, midplane, 1e-12));
    }
}

TEST_CASE("cutting_off_centre_splits_the_volume_in_the_same_ratio", "[BlockTestSuite]") {
    // At exactly 0.5 a swapped side, a transposed frame and a mirrored half all still measure 1 and
    // 1. Off-centre they cannot, which makes this the sharper of the 2 checks.
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_hex_block(tall_box());

    const auto e = edge_between(blocking, Point3d(0.0, 0.0, 0.0), Point3d(0.0, 0.0, 2.0));
    REQUIRE(blocking.cut_sheet(e, 0.25));

    auto volumes = blocking.block_volumes(4);
    REQUIRE(volumes.size() == 2);
    std::sort(volumes.begin(), volumes.end());
    REQUIRE(volumes[0] == Approx(0.5).margin(1e-12));
    REQUIRE(volumes[1] == Approx(1.5).margin(1e-12));

    // A quarter of the way from the edge's own start, so the cut plane sits at z=0.5 and every
    // corner of the blocking lies on one of the 3 planes it leaves behind.
    for (auto it = blocking.cmap().attributes<0>().begin(), end = blocking.cmap().attributes<0>().end(); it != end;
         ++it) {
        const double z = it->info().point.z();
        REQUIRE((std::abs(z) < 1e-12 || std::abs(z - 0.5) < 1e-12 || std::abs(z - 2.0) < 1e-12));
    }
}

TEST_CASE("every_edge_of_a_block_cuts_it_into_2_positive_volumes_summing_to_the_whole", "[BlockTestSuite]") {
    // Sweeps all 12 edges rather than trusting one: which frame each half lands in is CGAL's choice,
    // and the frame bugs this guards against showed on some edges while sparing others.
    for (int index = 0; index < 12; ++index) {
        const FacetedGeometry geom = make_far_geom_model();
        Blocking<FacetedGeometry> blocking(geom);
        blocking.create_hex_block(tall_box());

        int seen = 0;
        Blocking<FacetedGeometry>::Edge target{};
        for (auto it = blocking.cmap().attributes<1>().begin(), end = blocking.cmap().attributes<1>().end(); it != end;
             ++it, ++seen) {
            if (seen == index) target = it;
        }
        REQUIRE(blocking.cut_sheet(target, 0.25));

        auto volumes = blocking.block_volumes(4);
        REQUIRE(volumes.size() == 2);
        // Positive, so neither half came out inverted, and adding back up to what they came from.
        REQUIRE(volumes[0] > 0.0);
        REQUIRE(volumes[1] > 0.0);
        REQUIRE(volumes[0] + volumes[1] == Approx(2.0).margin(1e-12));
        std::sort(volumes.begin(), volumes.end());
        REQUIRE(volumes[0] == Approx(0.5).margin(1e-12));
        REQUIRE(volumes[1] == Approx(1.5).margin(1e-12));
    }
}
