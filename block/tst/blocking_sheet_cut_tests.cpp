#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iterator>
#include <limits>
#include <string>
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

    /** @brief The edge of @p ABlocking that runs along @p AAxis and whose midpoint is nearest
     * @p APoint. Positions are only ever compared for *nearness* here, never for equality: a cut
     * places its node by evaluating a curve, so the node it leaves at parameter 0.37 of a unit edge
     * is near 0.37 and not equal to it. */
    template<typename TBlocking>
    typename TBlocking::Edge edge_running_near(TBlocking &ABlocking, int AAxis, const Point3d &APoint) {
        auto &map = ABlocking.cmap();
        typename TBlocking::Edge best{};
        double best_d = 1e30;
        for (auto it = map.template attributes<1>().begin(), end = map.template attributes<1>().end(); it != end;
             ++it) {
            const auto d = it->dart();
            const Point3d &a = map.template attribute<0>(d)->info().point;
            const Point3d &b = map.template attribute<0>(map.template beta<1>(d))->info().point;
            const std::array<double, 3> pa{a.x(), a.y(), a.z()};
            const std::array<double, 3> pb{b.x(), b.y(), b.z()};
            bool runs_along = true;
            for (int k = 0; k < 3; ++k) {
                const bool same = std::abs(pa[k] - pb[k]) < 1e-9;
                if ((k == AAxis) == same) runs_along = false;
            }
            if (!runs_along) continue;
            const double dist = Vector3d(a + Vector3d(a, b) * 0.5, APoint).norm();
            if (dist < best_d) {
                best_d = dist;
                best = it;
            }
        }
        REQUIRE(best_d < 1e30);
        return best;
    }

    /** @brief How far the straightest possible edge would have to bend to become @p ACurve — the
     * largest distance from a sampled point on it to the chord joining its 2 ends. On a blocking cut
     * out of a box every edge is straight, so anything above rounding is a defect. */
    template<typename TCurve>
    double bend_of(const TCurve &ACurve) {
        const Point3d a = ACurve.value(0.0);
        const Point3d b = ACurve.value(1.0);
        const Vector3d chord(a, b);
        const double len2 = chord.dot(chord);
        double worst = 0.0;
        for (int i = 1; i < 16; ++i) {
            const Point3d p = ACurve.value(i / 16.0);
            const double t = (len2 > 0.0) ? std::clamp(Vector3d(a, p).dot(chord) / len2, 0.0, 1.0) : 0.0;
            worst = std::max(worst, Vector3d(p, a + chord * t).norm());
        }
        return worst;
    }

    /** @brief The worst `bend_of()` over every edge of @p ABlocking. */
    template<typename TBlocking>
    double worst_edge_bend(TBlocking &ABlocking) {
        double worst = 0.0;
        auto &map = ABlocking.cmap();
        for (auto it = map.template attributes<1>().begin(), end = map.template attributes<1>().end(); it != end;
             ++it) {
            worst = std::max(worst, bend_of(it->info().curve));
        }
        return worst;
    }

    /** @brief The block of @p ABlocking whose corner centroid is nearest @p APoint — a way to name a
     * block for deletion that does not depend on traversal order. */
    template<typename TBlocking>
    typename TBlocking::Block block_near(TBlocking &ABlocking, const Point3d &APoint) {
        auto &map = ABlocking.cmap();
        typename TBlocking::Block best{};
        double best_d = 1e30;
        for (auto it = map.template attributes<3>().begin(), end = map.template attributes<3>().end(); it != end;
             ++it) {
            Vector3d acc(0.0, 0.0, 0.0);
            int n = 0;
            const auto d0 = it->dart();
            const Point3d &origin = map.template attribute<0>(d0)->info().point;
            for (auto c = map.template one_dart_per_incident_cell<0, 3>(d0).begin(),
                      cend = map.template one_dart_per_incident_cell<0, 3>(d0).end();
                 c != cend;
                 ++c, ++n) {
                acc += Vector3d(origin, map.template attribute<0>(c)->info().point);
            }
            const Point3d centre = origin + acc * (1.0 / static_cast<double>(n));
            const double d = Vector3d(centre, APoint).norm();
            if (d < best_d) {
                best_d = d;
                best = it;
            }
        }
        return best;
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
    // Curved the way a user gets there — a block fitted onto a real cylinder — rather than by
    // bending control points by hand. That matters: classify() leaves each cell's curve, surface and
    // volume agreeing with one another, and a cut relies on exactly that agreement when it pins a
    // new edge onto the nodes the sheet placed.
    std::string dir(TEST_SAMPLES_DIR);
    const FacetedGeometry geom(dir + "/cylinder.msh");
    using CubicBlocking = Blocking<FacetedGeometry>;

    const auto fitted = [&geom](CubicBlocking &ABlocking) {
        // A box around the model, then snapped onto it: the corners land on the cylinder and the
        // edges bend to follow it.
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
        ABlocking.create_hex_block({Point3d(lo[0], lo[1], lo[2]),
                                    Point3d(hi[0], lo[1], lo[2]),
                                    Point3d(hi[0], hi[1], lo[2]),
                                    Point3d(lo[0], hi[1], lo[2]),
                                    Point3d(lo[0], lo[1], hi[2]),
                                    Point3d(hi[0], lo[1], hi[2]),
                                    Point3d(hi[0], hi[1], hi[2]),
                                    Point3d(lo[0], hi[1], hi[2])});
        ABlocking.classify(0.3);
    };

    CubicBlocking reference(geom);
    fitted(reference);
    CubicBlocking cut(geom);
    fitted(cut);

    // Guard against silently testing a straight block: fitting has to reach the sampled geometry.
    const auto bent_mesh = reference.to_mesh(3);
    double largest_shift = 0.0;
    {
        CubicBlocking straight(geom);
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
        straight.create_hex_block({Point3d(lo[0], lo[1], lo[2]),
                                   Point3d(hi[0], lo[1], lo[2]),
                                   Point3d(hi[0], hi[1], lo[2]),
                                   Point3d(lo[0], hi[1], lo[2]),
                                   Point3d(lo[0], lo[1], hi[2]),
                                   Point3d(hi[0], lo[1], hi[2]),
                                   Point3d(hi[0], hi[1], hi[2]),
                                   Point3d(lo[0], hi[1], hi[2])});
        const auto straight_mesh = straight.to_mesh(3);
        REQUIRE(bent_mesh.nb_nodes() == straight_mesh.nb_nodes());
        for (UInt i = 0; i < bent_mesh.nb_nodes(); ++i) {
            largest_shift =
                std::max(largest_shift, Vector3d(bent_mesh.node(NodeId{i}), straight_mesh.node(NodeId{i})).norm());
        }
    }
    REQUIRE(largest_shift > 1e-2);

    const auto fine = reference.to_mesh(12);
    std::vector<Point3d> cloud;
    cloud.reserve(fine.nb_nodes());
    for (UInt i = 0; i < fine.nb_nodes(); ++i) {
        cloud.push_back(fine.node(NodeId{i}));
    }

    // Cut at 1/3 and mesh each half with 4 intervals: that samples the parent at 0,1/12,...,12/12
    // along the cut axis and at 0,3/12,...,12/12 across the other 2 — every one a node the uncut
    // blocking's own 12-interval mesh already carries. So the cut mesh must be a subset of it.
    int seen = 0;
    CubicBlocking::Edge target{};
    for (auto it = cut.cmap().attributes<1>().begin(), end = cut.cmap().attributes<1>().end(); it != end;
         ++it, ++seen) {
        if (seen == 0) target = it;
    }
    REQUIRE(cut.cut_sheet(target, 1.0 / 3.0));
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

TEST_CASE("a_sheet_stays_in_one_plane_after_an_earlier_cut", "[BlockTestSuite]") {
    // The failure this guards against needed 3 things at once, which is why every earlier test in
    // this file missed it: a *second* cut (the first has no cut-made edges yet), an *off-centre*
    // parameter (at 0.5, cutting an edge from the wrong end lands in the same place), and an edge
    // the previous cut created inside a face — whose curve comes from a surface's control grid while
    // its nodes were placed by evaluating another curve. Those agree as real numbers and not as
    // doubles, so deciding which end a curve starts at by exact position equality answers "the other
    // one" for exactly those edges, cutting them at 1-t and pulling the sheet out of plane.
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_hex_block(tall_box());

    const auto e = edge_between(blocking, Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0));
    REQUIRE(blocking.cut_sheet(e, 0.4));

    // Every edge of the once-cut blocking, not just the ones the next cut happens to aim at.
    int checked = 0;
    for (auto it = blocking.cmap().attributes<1>().begin(), end = blocking.cmap().attributes<1>().end(); it != end;
         ++it) {
        const auto sheet = blocking.find_sheet(it);
        REQUIRE(sheet.has_value());

        // The cut runs perpendicular to the edges it splits, so on this axis-aligned block every
        // cut point must share the coordinate of the axis the aimed edge runs along.
        const auto d = it->dart();
        const Point3d &a = blocking.cmap().attribute<0>(d)->info().point;
        const Point3d &b = blocking.cmap().attribute<0>(blocking.cmap().beta<1>(d))->info().point;
        const std::array<double, 3> pa{a.x(), a.y(), a.z()};
        const std::array<double, 3> pb{b.x(), b.y(), b.z()};
        std::size_t axis = 0;
        for (std::size_t k = 0; k < 3; ++k) {
            if (std::abs(pa[k] - pb[k]) > std::abs(pa[axis] - pb[axis])) axis = k;
        }

        double lo = std::numeric_limits<double>::max();
        double hi = std::numeric_limits<double>::lowest();
        for (const auto &se : sheet->edges) {
            const Point3d p = blocking.cut_point(se, 0.4);
            const std::array<double, 3> c{p.x(), p.y(), p.z()};
            lo = std::min(lo, c[axis]);
            hi = std::max(hi, c[axis]);
        }
        REQUIRE(hi - lo == Approx(0.0).margin(1e-12));
        ++checked;
    }
    REQUIRE(checked == 20);
}

TEST_CASE("repeated_off_centre_cuts_keep_the_generated_mesh_on_the_original_block", "[BlockTestSuite]") {
    // The same failure seen from where a user would see it: the mesh. A sheet cut out of plane still
    // leaves a valid map with sensible block volumes, and only the generated mesh — the thing biy
    // actually draws — shows it, by drifting off the block it came from.
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_hex_block(tall_box());

    const auto measure = [&]() {
        const auto mesh = blocking.to_mesh(2);
        double lo_x = 0.0;
        double hi_x = 0.0;
        double hi_z = 0.0;
        for (UInt i = 0; i < mesh.nb_nodes(); ++i) {
            const Point3d &p = mesh.node(NodeId{i});
            lo_x = std::min(lo_x, p.x());
            hi_x = std::max(hi_x, p.x());
            hi_z = std::max(hi_z, p.z());
        }
        // Nothing the cut generates may leave the box the blocking started as.
        REQUIRE(lo_x >= Approx(0.0).margin(1e-9));
        REQUIRE(hi_x <= Approx(1.0).margin(1e-9));
        REQUIRE(hi_z <= Approx(2.0).margin(1e-9));
    };

    measure();
    for (const auto &pair : {std::pair{Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0)},
                             std::pair{Point3d(0.0, 0.0, 0.0), Point3d(0.0, 1.0, 0.0)},
                             std::pair{Point3d(0.0, 0.0, 0.0), Point3d(0.0, 0.0, 2.0)}}) {
        const auto target = edge_between(blocking, pair.first, pair.second);
        REQUIRE(blocking.cut_sheet(target, 0.4));
        REQUIRE(blocking.is_valid_topology());
        measure();

        // And the blocking still adds back up to what it started as, block by block.
        const auto volumes = blocking.block_volumes(2);
        double total = 0.0;
        for (const double v : volumes) {
            REQUIRE(v > 0.0);
            total += v;
        }
        REQUIRE(total == Approx(2.0).margin(1e-12));
    }
    REQUIRE(blocking.template nb_cells<3>() == 8);
}

TEST_CASE("deleting_a_lone_block_leaves_nothing_behind", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_hex_block(box(0.0, 1.0));
    REQUIRE(blocking.template nb_cells<3>() == 1);

    blocking.delete_block(blocking.cmap().attributes<3>().begin());

    // Everything the block owned went with it: nothing was shared, so nothing survives.
    REQUIRE(blocking.template nb_cells<3>() == 0);
    REQUIRE(blocking.template nb_cells<2>() == 0);
    REQUIRE(blocking.template nb_cells<1>() == 0);
    REQUIRE(blocking.template nb_cells<0>() == 0);
    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.to_mesh(1).nb_nodes() == 0);
}

TEST_CASE("deleting_one_of_2_sewn_blocks_keeps_what_they_shared", "[BlockTestSuite]") {
    // The case that says whether the removal cascades correctly: the shared face, its 4 edges and
    // its 4 corners must survive as the survivor's own boundary, while the 5 faces, 8 edges and 4
    // corners that belonged to the deleted block alone must go.
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_hex_block(box(0.0, 1.0));
    blocking.create_hex_block(box(1.0, 2.0));
    blocking.build_connectivity();
    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.template nb_cells<0>() == 12);
    REQUIRE(blocking.template nb_cells<2>() == 11);

    blocking.delete_block(blocking.cmap().attributes<3>().begin());

    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.template nb_cells<3>() == 1);
    REQUIRE(blocking.template nb_cells<0>() == 8); // one cube's worth
    REQUIRE(blocking.template nb_cells<1>() == 12);
    REQUIRE(blocking.template nb_cells<2>() == 6);

    // And what is left still meshes as one whole hex, the shared face now being a boundary one.
    const auto mesh = blocking.to_mesh(1);
    REQUIRE(mesh.nb_cells() == 1);
    REQUIRE(mesh.nb_nodes() == 8);
    REQUIRE(mesh.nb_faces() == 0); // a block's faces are not emitted as standalone quads
}

TEST_CASE("deleting_a_block_after_a_cut_keeps_the_face_index_sound", "[BlockTestSuite]") {
    // A cut leaves faces belonging to blocks, and a deletion garbage-collects some of them. Once
    // CGAL has collected a face attribute there is no handle left to look up and drop, so the
    // "which faces belong to a block" index has to be re-derived — if it were not, to_mesh() would
    // compare against freed attributes and emit stray quads for the survivor's own faces.
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_hex_block(box(0.0, 1.0));

    const auto e = edge_along(blocking, 0, 0.0);
    REQUIRE(blocking.cut_sheet(e, 0.4));
    REQUIRE(blocking.template nb_cells<3>() == 2);

    blocking.delete_block(blocking.cmap().attributes<3>().begin());
    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.template nb_cells<3>() == 1);

    const auto mesh = blocking.to_mesh(2);
    REQUIRE(mesh.nb_cells() == 8);
    // No standalone quads: every remaining face belongs to the surviving block.
    REQUIRE(mesh.nb_faces() == 0);

    // The survivor is still a genuine block, and can be cut again.
    const auto again = blocking.cmap().attributes<1>().begin();
    REQUIRE(blocking.find_sheet(again).has_value());
}

TEST_CASE("deleting_every_block_one_by_one_empties_the_blocking", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    blocking.create_hex_block(box(0.0, 1.0));
    blocking.create_hex_block(box(1.0, 2.0));
    blocking.build_connectivity();

    while (blocking.template nb_cells<3>() > 0) {
        blocking.delete_block(blocking.cmap().attributes<3>().begin());
        REQUIRE(blocking.is_valid_topology());
    }
    REQUIRE(blocking.template nb_cells<2>() == 0);
    REQUIRE(blocking.template nb_cells<0>() == 0);

    // An empty blocking is a legitimate state to be in, not a broken one: it still builds a mesh
    // and still takes a new block.
    REQUIRE(blocking.to_mesh(2).nb_nodes() == 0);
    blocking.create_hex_block(box(0.0, 1.0));
    REQUIRE(blocking.template nb_cells<3>() == 1);
    REQUIRE(blocking.to_mesh(1).nb_cells() == 1);
}

TEST_CASE("cutting_after_a_deletion_still_cuts_every_edge_from_the_same_end", "[BlockTestSuite]") {
    // Cut a box, take 2 blocks out of it, and cut again. Every edge is straight throughout — it is
    // an axis-aligned box cut by axis-aligned sheets — so any bend at all means a node was placed
    // somewhere it does not belong.
    //
    // What used to put one there: an edge names the end its stored curve starts from, and it named
    // it by comparing the 2 endpoint attribute *handles*, which are addresses. Deleting a block has
    // CGAL rebuild the attributes whose vertex orbit came apart, and the 2 addresses can then come
    // back in the opposite order — while the stored curve, of course, has not moved. The edge then
    // reports the end it does not start at, so the next sheet cut through it lands at `1-t` there
    // and at `t` on all its parallels, leaving one node out of line and every edge through it bent.
    //
    // Cells are named by their position in the map's own traversal, and the sequence is the shortest
    // one a search over random cut/delete sequences found that trips it: whether a given deletion
    // flips a given pair of addresses is not something a geometric description of the operations can
    // pin down, so a readable "cut this edge, delete that block" sequence does not reproduce it.
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_hex_block(box(0.0, 1.0));

    const auto cut_nth = [&](std::size_t AN, double AT) {
        auto it = blocking.cmap().attributes<1>().begin();
        std::advance(it, static_cast<long>(AN));
        REQUIRE(blocking.cut_sheet(it, AT));
        REQUIRE(blocking.is_valid_topology());
        REQUIRE(worst_edge_bend(blocking) < 1e-12);
    };
    const auto delete_nth = [&](std::size_t AN) {
        auto it = blocking.cmap().attributes<3>().begin();
        std::advance(it, static_cast<long>(AN));
        blocking.delete_block(it);
        REQUIRE(blocking.is_valid_topology());
        REQUIRE(worst_edge_bend(blocking) < 1e-12);
    };

    cut_nth(10, 0.6);
    cut_nth(3, 0.221);
    delete_nth(3);
    delete_nth(2);
    cut_nth(20, 0.239);
}

TEST_CASE("cutting_a_holed_grid_leaves_every_edge_straight", "[BlockTestSuite]") {
    // The same property stated the way a user meets it, and independent of any one traversal order:
    // cut a box into a grid, punch holes in it, keep cutting, and nothing ever bends.
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_hex_block(box(0.0, 1.0));

    const auto cut_along = [&](int AAxis, const Point3d &ANear, double AT) {
        const auto target = edge_running_near(blocking, AAxis, ANear);
        REQUIRE(blocking.cut_sheet(target, AT));
        REQUIRE(blocking.is_valid_topology());
        REQUIRE(worst_edge_bend(blocking) < 1e-12);
    };
    const auto delete_near = [&](const Point3d &ANear) {
        blocking.delete_block(block_near(blocking, ANear));
        REQUIRE(blocking.is_valid_topology());
        REQUIRE(worst_edge_bend(blocking) < 1e-12);
    };

    // A 2x2x2 grid, cut off-centre so a mirrored cut cannot pass unnoticed.
    cut_along(0, Point3d(0.5, 0.0, 0.0), 0.37);
    cut_along(1, Point3d(0.0, 0.5, 0.0), 0.43);
    cut_along(2, Point3d(0.0, 0.0, 0.5), 0.61);
    REQUIRE(blocking.template nb_cells<3>() == 8);

    delete_near(Point3d(0.1, 0.1, 0.1));
    delete_near(Point3d(0.9, 0.9, 0.9));
    REQUIRE(blocking.template nb_cells<3>() == 6);

    cut_along(0, Point3d(0.7, 0.0, 0.0), 0.29);
    cut_along(1, Point3d(0.0, 0.7, 0.0), 0.71);
    cut_along(2, Point3d(0.0, 0.0, 0.8), 0.53);

    delete_near(Point3d(0.9, 0.1, 0.1));
    delete_near(Point3d(0.1, 0.9, 0.9));

    cut_along(0, Point3d(0.2, 0.0, 0.0), 0.44);
    cut_along(1, Point3d(0.0, 0.2, 0.0), 0.36);
    cut_along(2, Point3d(0.0, 0.0, 0.3), 0.62);
}

TEST_CASE("deleting_a_sheet_glues_back_what_was_either_side_of_it", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_hex_block(box(0.0, 1.0));

    // 3 slabs along z, then the middle one taken out.
    REQUIRE(blocking.cut_sheet(edge_running_near(blocking, 2, Point3d(0.0, 0.0, 0.5)), 1.0 / 3.0));
    REQUIRE(blocking.cut_sheet(edge_running_near(blocking, 2, Point3d(0.0, 0.0, 0.66)), 0.5));
    REQUIRE(blocking.template nb_cells<3>() == 3);

    const auto middle = edge_running_near(blocking, 2, Point3d(0.0, 0.0, 0.5));
    REQUIRE(blocking.delete_sheet(middle, 1e-9));

    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.template nb_cells<3>() == 2);
    REQUIRE(worst_edge_bend(blocking) < 1e-12);

    // The 2 survivors still fill the box they started as.
    double total = 0.0;
    for (const double v : blocking.block_volumes(2)) {
        REQUIRE(v > 0.0);
        total += v;
    }
    REQUIRE(total == Approx(1.0).margin(1e-12));
}

TEST_CASE("deleting_a_boundary_sheet_pulls_that_side_of_the_blocking_in", "[BlockTestSuite]") {
    // A sheet with blocks on only one side of it: collapsing it has nothing to glue on the far side,
    // so that face of the blocking simply moves in to where the merged corners land.
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_hex_block(box(0.0, 1.0));
    REQUIRE(blocking.cut_sheet(edge_running_near(blocking, 0, Point3d(0.5, 0.0, 0.0)), 0.5));
    REQUIRE(blocking.cut_sheet(edge_running_near(blocking, 1, Point3d(0.0, 0.5, 0.0)), 0.5));
    REQUIRE(blocking.cut_sheet(edge_running_near(blocking, 2, Point3d(0.0, 0.0, 0.5)), 0.5));
    REQUIRE(blocking.template nb_cells<3>() == 8);

    // The bottom layer, whose corners at z=0 and z=0.5 merge at z=0.25.
    REQUIRE(blocking.delete_sheet(edge_running_near(blocking, 2, Point3d(0.0, 0.0, 0.25)), 1e-9));

    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.template nb_cells<3>() == 4);
    REQUIRE(worst_edge_bend(blocking) < 1e-12);
    double total = 0.0;
    for (const double v : blocking.block_volumes(2)) {
        REQUIRE(v > 0.0);
        total += v;
    }
    // What is left spans z in [0.25, 1].
    REQUIRE(total == Approx(0.75).margin(1e-12));
}

TEST_CASE("deleting_a_sheet_contracts_the_blocking_rather_than_undoing_the_cut", "[BlockTestSuite]") {
    // Worth being explicit about, because it is the first thing one expects and it is not what
    // happens: collapsing the layer a cut just created does *not* put the block back. The cut adds a
    // layer between 2 existing faces; collapsing one pulls its 2 faces together onto the merge, so a
    // box cut at 0.4 of its height and then collapsed on the lower layer keeps only what is above
    // the merged corners. 1 block again, but a shorter one.
    //
    // On a classified blocking the merge lands differently: a corner on the model wins over an
    // interior one and is projected back onto it, so a boundary stays where the model puts it (see
    // the classification tests). Here nothing is classified, so the merge is simply the midpoint.
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_hex_block(tall_box());
    REQUIRE(blocking.cut_sheet(edge_running_near(blocking, 2, Point3d(0.0, 0.0, 1.0)), 0.4));
    REQUIRE(blocking.template nb_cells<3>() == 2);

    REQUIRE(blocking.delete_sheet(edge_running_near(blocking, 2, Point3d(0.0, 0.0, 0.4)), 1e-9));

    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.template nb_cells<3>() == 1);
    REQUIRE(blocking.template nb_cells<2>() == 6);
    REQUIRE(blocking.template nb_cells<1>() == 12);
    REQUIRE(blocking.template nb_cells<0>() == 8);
    REQUIRE(worst_edge_bend(blocking) < 1e-12);
    // z=0 and z=0.8 merged at z=0.4, leaving the box from there to its top at z=2.
    REQUIRE(blocking.block_volumes(2)[0] == Approx(1.6).margin(1e-12));
    double lowest = 1e30;
    auto &map = blocking.cmap();
    for (auto it = map.attributes<0>().begin(), end = map.attributes<0>().end(); it != end; ++it) {
        lowest = std::min(lowest, it->info().point.z());
    }
    REQUIRE(lowest == Approx(0.4).margin(1e-12));
}

TEST_CASE("deleting_a_sheet_leaves_the_blocks_it_ran_between_sharing_a_face", "[BlockTestSuite]") {
    // The gluing, seen on 2 blocks that were sewn along the collapsed layer's far side: after it
    // goes, they are 1 block, and the boundary it used to sit against has moved onto the merge.
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_hex_block(box(0.0, 1.0));
    blocking.create_hex_block(box(1.0, 2.0));
    blocking.build_connectivity();
    REQUIRE(blocking.template nb_cells<3>() == 2);

    // The sheet through an x-running edge of the first block is that block's own 4: it does not
    // reach the second, whose x-edges are separated from them by the face the 2 share.
    REQUIRE(blocking.delete_sheet(edge_running_near(blocking, 0, Point3d(0.5, 0.0, 0.0)), 1e-9));

    REQUIRE(blocking.is_valid_topology());
    REQUIRE(blocking.template nb_cells<3>() == 1);
    REQUIRE(blocking.template nb_cells<0>() == 8);
    REQUIRE(worst_edge_bend(blocking) < 1e-12);
    // The survivor now runs from the merged corners at x=0.5 to x=2.
    REQUIRE(blocking.block_volumes(2)[0] == Approx(1.5).margin(1e-12));
}

TEST_CASE("delete_sheet_refuses_a_sheet_that_is_the_whole_blocking", "[BlockTestSuite]") {
    // Nothing either side to glue: collapsing would empty the blocking rather than thin it, so it is
    // refused and nothing moves.
    const FacetedGeometry geom = make_far_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 3);
    blocking.create_hex_block(box(0.0, 1.0));

    REQUIRE_FALSE(blocking.delete_sheet(edge_running_near(blocking, 2, Point3d(0.0, 0.0, 0.5)), 1e-9));

    REQUIRE(blocking.template nb_cells<3>() == 1);
    REQUIRE(blocking.template nb_cells<2>() == 6);
    REQUIRE(blocking.template nb_cells<1>() == 12);
    REQUIRE(blocking.template nb_cells<0>() == 8);
    REQUIRE(blocking.block_volumes(2)[0] == Approx(1.0).margin(1e-12));

    // And on 2 blocks the sheet crosses both, which is again all of them.
    blocking.create_hex_block(box(1.0, 2.0));
    blocking.build_connectivity();
    REQUIRE_FALSE(blocking.delete_sheet(edge_running_near(blocking, 2, Point3d(0.0, 0.0, 0.5)), 1e-9));
    REQUIRE(blocking.template nb_cells<3>() == 2);
}
