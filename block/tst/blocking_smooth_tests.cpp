#include <algorithm>
#include <array>
#include <filesystem>
#include <vector>

#include <gecko/block/Blocking.h>
#include <gecko/block/Smoother.h>
#include <gecko/geom/FacetedGeometry.h>
#include <gecko/io/GmshMeshWriter.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
using Catch::Approx;

using namespace gecko;

namespace {
    /**
     * @brief The unit square as a geometric model: 4 tagged vertices (1..4), 4 tagged boundary
     * curves (10..13, bottom/right/top/left) and one tagged surface (20).
     *
     * The same fixture `blocking_classification_tests.cpp` uses, and for the same reason: it is the
     * smallest model with something at every dimension, so a blocking laid on it has corners pinned
     * to vertices, corners riding curves and corners free on a surface all at once.
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
        const std::array<NodeId, 4> vertices{v0, v1, v2, v3};
        for (std::size_t v = 0; v < 4; ++v) {
            node_group[vertices[v].value] = vtx_group;
            node_entity[vertices[v].value] = static_cast<Int>(v) + 1;
        }

        auto &edge_group = mesh.add_variable<GroupId, CellType::Edge>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &edge_entity = mesh.add_variable<Int, CellType::Edge>(std::string(io::ENTITY_TAG_VARIABLE));
        const std::array<std::pair<NodeId, NodeId>, 4> sides{
            std::pair{v0, v1}, std::pair{v1, v2}, std::pair{v2, v3}, std::pair{v3, v0}};
        for (std::size_t s = 0; s < 4; ++s) {
            auto e = mesh.add_edge(sides[s].first, sides[s].second);
            edge_group[e.value] = curve_group;
            edge_entity[e.value] = 10 + static_cast<Int>(s);
        }

        auto &face_group = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &face_entity = mesh.add_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
        for (const auto &tri : {std::array{v0, v1, v2}, std::array{v0, v2, v3}}) {
            auto f = mesh.add_face(tri[0], tri[1], tri[2]);
            face_group[f.value] = surf_group;
            face_entity[f.value] = 20;
        }

        const auto path = (std::filesystem::temp_directory_path() / "gecko_block_smooth_square_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }

    /** @brief Builds the 2x2 grid of quad blocks tiling the unit square, sewn together. */
    void build_2x2_grid(Blocking<FacetedGeometry> &ABlocking) {
        for (const double x : {0.0, 0.5}) {
            for (const double y : {0.0, 0.5}) {
                ABlocking.create_quad_block({Point3d(x, y, 0.0),
                                             Point3d(x + 0.5, y, 0.0),
                                             Point3d(x + 0.5, y + 0.5, 0.0),
                                             Point3d(x, y + 0.5, 0.0)});
            }
        }
        ABlocking.build_connectivity();
    }

    /** @brief Builds the 2x2x2 grid of unit hex blocks filling [0,2]^3, sewn together. */
    void build_2x2x2_grid(Blocking<FacetedGeometry> &ABlocking) {
        for (const double x : {0.0, 1.0}) {
            for (const double y : {0.0, 1.0}) {
                for (const double z : {0.0, 1.0}) {
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

    /** @brief Finds the node sitting exactly at @p APoint. */
    auto node_at(Blocking<FacetedGeometry> &ABlocking, const Point3d &APoint) {
        auto &map = ABlocking.cmap();
        for (auto it = map.attributes<0>().begin(), itend = map.attributes<0>().end(); it != itend; ++it) {
            if (it->info().point == APoint) return it;
        }
        FAIL("no node found at the requested position");
        return map.attributes<0>().begin();
    }

    /** @brief Writes a corner's position straight into the map, the way a user's edit leaves it —
     * without move_node(), so nothing around it is rebuilt and the smoother has to do it. */
    void displace(Blocking<FacetedGeometry> &ABlocking, const Point3d &AFrom, const Point3d &ATo) {
        node_at(ABlocking, AFrom)->info().point = ATo;
    }

    /** @brief How far a blocking's corners sit, in total, from a recorded set of positions —
     * compared as sorted sets, since nothing promises the traversal order matches. */
    double distance_from(Blocking<FacetedGeometry> &ABlocking, const std::vector<std::array<double, 3>> &AReference);

    /** @brief Every corner position, sorted, so 2 blockings can be compared as sets of corners. */
    std::vector<std::array<double, 3>> sorted_positions(Blocking<FacetedGeometry> &ABlocking) {
        std::vector<std::array<double, 3>> points;
        auto &map = ABlocking.cmap();
        for (auto it = map.attributes<0>().begin(), itend = map.attributes<0>().end(); it != itend; ++it) {
            points.push_back({it->info().point.x(), it->info().point.y(), it->info().point.z()});
        }
        std::ranges::sort(points);
        return points;
    }

    double distance_from(Blocking<FacetedGeometry> &ABlocking, const std::vector<std::array<double, 3>> &AReference) {
        const auto points = sorted_positions(ABlocking);
        double total = 0.0;
        for (std::size_t p = 0; p < points.size() && p < AReference.size(); ++p) {
            for (std::size_t c = 0; c < 3; ++c) {
                const double d = points[p][c] - AReference[p][c];
                total += d * d;
            }
        }
        return total;
    }
} // namespace

using Smooth = Smoother<FacetedGeometry>;

TEST_CASE("smoothing_pulls_a_displaced_interior_corner_back_to_the_average_of_its_neighbours", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    build_2x2_grid(blocking);
    blocking.classify(0.05);
    displace(blocking, Point3d(0.5, 0.5, 0.0), Point3d(0.35, 0.6, 0.0));

    Smooth smoother(blocking);
    const auto report = smoother.smooth(10);

    // Only the centre had anywhere to go: every other corner is on a model vertex (pinned) or is
    // already at the midpoint of its own curve.
    REQUIRE(report.moves == 1);
    const auto centre = node_at(blocking, Point3d(0.5, 0.5, 0.0));
    REQUIRE(centre->info().point.x() == Approx(0.5));
    REQUIRE(centre->info().point.y() == Approx(0.5));
    REQUIRE(centre->info().point.z() == Approx(0.0).margin(1e-15));
}

TEST_CASE("smoothing_puts_a_whole_perturbed_grid_back_the_way_it_was", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    build_2x2_grid(blocking);
    blocking.classify(0.05);
    const auto regular = sorted_positions(blocking);

    displace(blocking, Point3d(0.5, 0.5, 0.0), Point3d(0.3, 0.65, 0.0));
    displace(blocking, Point3d(0.5, 0.0, 0.0), Point3d(0.2, 0.0, 0.0));
    displace(blocking, Point3d(0.0, 0.5, 0.0), Point3d(0.0, 0.8, 0.0));
    REQUIRE(distance_from(blocking, regular) > 0.0);

    Smooth smoother(blocking);
    const auto report = smoother.smooth(50);

    REQUIRE(distance_from(blocking, regular) == Approx(0.0).margin(1e-18));
    REQUIRE(report.worst_quality == Approx(1.0));
}

TEST_CASE("the_optimization_pass_improves_on_where_the_laplacian_stops", "[BlockTestSuite]") {
    // A corner of the boundary held off-centre, so the average of the 4 neighbours is no longer
    // where the 4 quads are happiest. The Laplacian has only that one position to offer and stops
    // there; the search carries on from it.
    const auto run = [](Smooth::Strategy AStrategy) {
        const FacetedGeometry geom = make_square_geom_model();
        Blocking<FacetedGeometry> blocking(geom);
        build_2x2_grid(blocking);
        blocking.classify(0.05);
        displace(blocking, Point3d(0.5, 0.0, 0.0), Point3d(0.15, 0.0, 0.0));

        Smooth smoother(blocking);
        smoother.lock(node_at(blocking, Point3d(0.15, 0.0, 0.0))->info().id);
        return smoother.smooth(60, AStrategy);
    };

    const auto laplacian = run(Smooth::Strategy::Laplacian);
    const auto both = run(Smooth::Strategy::Both);

    REQUIRE(laplacian.optimization_passes == 0);
    REQUIRE(both.laplacian_passes == laplacian.laplacian_passes);
    REQUIRE(both.optimization_passes > 0);
    REQUIRE(both.worst_quality > laplacian.worst_quality);
}

TEST_CASE("smoothing_leaves_a_corner_pinned_to_a_model_vertex_alone", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    build_2x2_grid(blocking);
    blocking.classify(0.05);
    displace(blocking, Point3d(0.5, 0.5, 0.0), Point3d(0.2, 0.75, 0.0));

    Smooth smoother(blocking);
    smoother.smooth(20);

    for (const Point3d &corner :
         {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)}) {
        const auto node = node_at(blocking, corner);
        REQUIRE(node->info().geom_targets.size() == 1);
        REQUIRE(node->info().geom_targets[0].first == GroupDim::Dim0);
    }
}

TEST_CASE("smoothing_keeps_a_corner_on_the_curve_it_is_classified_on", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    build_2x2_grid(blocking);
    blocking.classify(0.05);
    displace(blocking, Point3d(0.5, 0.0, 0.0), Point3d(0.15, 0.0, 0.0));

    Smooth smoother(blocking);
    smoother.smooth(20);

    // Slid back along the bottom curve towards its middle, and never off it: the interior corner
    // above it is not one of the neighbours it is averaged over, so nothing pulls it up into the
    // surface in the first place.
    const auto node = node_at(blocking, Point3d(0.5, 0.0, 0.0));
    REQUIRE(node->info().point.x() == Approx(0.5));
    REQUIRE(node->info().point.y() == Approx(0.0).margin(1e-15));
    REQUIRE(node->info().geom_targets == std::vector{std::pair{GroupDim::Dim1, Int(10)}});
}

TEST_CASE("smoothing_never_moves_a_locked_corner", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    build_2x2_grid(blocking);
    blocking.classify(0.05);
    displace(blocking, Point3d(0.5, 0.5, 0.0), Point3d(0.35, 0.6, 0.0));

    Smooth smoother(blocking);
    smoother.lock(node_at(blocking, Point3d(0.35, 0.6, 0.0))->info().id);
    const auto report = smoother.smooth(20);

    // Everything else is free to rearrange itself around it — and does, the corner it has to work
    // around being off. The one assertion that matters is that the locked corner is not among them.
    const auto centre = node_at(blocking, Point3d(0.35, 0.6, 0.0));
    REQUIRE(centre->info().point == Point3d(0.35, 0.6, 0.0));
}

TEST_CASE("smoothing_a_regular_grid_changes_nothing_and_stops_after_one_pass", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    build_2x2_grid(blocking);
    blocking.classify(0.05);
    const auto before = sorted_positions(blocking);

    Smooth smoother(blocking);
    const auto report = smoother.smooth(20);

    REQUIRE(report.moves == 0);
    REQUIRE(report.laplacian_passes == 1);
    REQUIRE(report.optimization_passes == 1);
    REQUIRE(sorted_positions(blocking) == before);
}

TEST_CASE("smoothing_never_makes_the_worst_cell_worse", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    build_2x2_grid(blocking);
    blocking.classify(0.05);
    displace(blocking, Point3d(0.5, 0.5, 0.0), Point3d(0.06, 0.94, 0.0));

    Smooth smoother(blocking);
    const double before = smoother.worst_quality();
    const auto report = smoother.smooth(30);

    // The guarantee that separates the smart Laplacian from the plain one, which would happily drag
    // this corner past its neighbours and turn a quad inside out.
    REQUIRE(report.worst_quality >= before);
    REQUIRE(report.worst_quality > 0.0);
}

TEST_CASE("smoothing_recomputes_the_control_points_of_every_cell_it_moved", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom, 3);
    build_2x2_grid(blocking);
    blocking.classify(0.05);
    displace(blocking, Point3d(0.5, 0.5, 0.0), Point3d(0.35, 0.6, 0.0));

    Smooth smoother(blocking);
    smoother.smooth(20);

    // The corner was written straight into the map, so before the refit every curve through it
    // still ran to where it used to be. Each edge's own curve has to end on its 2 corners again.
    auto &map = blocking.cmap();
    for (auto it = map.attributes<1>().begin(), itend = map.attributes<1>().end(); it != itend; ++it) {
        const auto dart = it->dart();
        const Point3d p0 = map.attribute<0>(dart)->info().point;
        const Point3d p1 = map.attribute<0>(map.beta<1>(dart))->info().point;
        const Point3d start = it->info().curve.value(0.0);
        const Point3d end = it->info().curve.value(1.0);
        const bool ends_on_its_corners = (start == p0 && end == p1) || (start == p1 && end == p0);
        REQUIRE(ends_on_its_corners);
    }
}

TEST_CASE("smoothing_works_the_same_on_a_3d_blocking", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_square_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    build_2x2x2_grid(blocking);

    // Unclassified, so nothing holds the outside in place: locked by hand instead, which is what
    // biy's frozen corners come to. Only the middle of the 27 is left free.
    Smooth smoother(blocking);
    auto &map = blocking.cmap();
    for (auto it = map.attributes<0>().begin(), itend = map.attributes<0>().end(); it != itend; ++it) {
        if (it->info().point != Point3d(1.0, 1.0, 1.0)) smoother.lock(it->info().id);
    }
    displace(blocking, Point3d(1.0, 1.0, 1.0), Point3d(1.45, 0.7, 1.3));

    const auto report = smoother.smooth(20);

    REQUIRE(report.moves >= 1);
    const auto centre = node_at(blocking, Point3d(1.0, 1.0, 1.0));
    REQUIRE(centre->info().point.x() == Approx(1.0));
    REQUIRE(centre->info().point.y() == Approx(1.0));
    REQUIRE(centre->info().point.z() == Approx(1.0));
    REQUIRE(report.worst_quality == Approx(1.0));
}
