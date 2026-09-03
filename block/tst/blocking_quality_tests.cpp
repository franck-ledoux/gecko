#include <array>
#include <cmath>
#include <filesystem>

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

        const auto path = (std::filesystem::temp_directory_path() / "gecko_block_quality_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }

    /** @brief The unit cube's 8 corners, in `HEX_CORNER_UVW` order. */
    std::array<Point3d, 8> unit_cube() {
        return {Point3d(0.0, 0.0, 0.0),
                Point3d(1.0, 0.0, 0.0),
                Point3d(1.0, 1.0, 0.0),
                Point3d(0.0, 1.0, 0.0),
                Point3d(0.0, 0.0, 1.0),
                Point3d(1.0, 0.0, 1.0),
                Point3d(1.0, 1.0, 1.0),
                Point3d(0.0, 1.0, 1.0)};
    }

    /** @brief The unit square's 4 corners in the z=0 plane, in `QUAD_CORNER_IJ` order. */
    std::array<Point3d, 4> unit_square() {
        return {Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)};
    }

    /** @brief Scales @p APoint by @p AScale, rotates it by @p AAngle about each axis in turn, then
     * translates it — an arbitrary similarity, which a shape measure must be blind to. */
    Point3d transformed(const Point3d &APoint, double AScale, double AAngle, const Vector3d &AShift) {
        const double c = std::cos(AAngle), s = std::sin(AAngle);
        double x = APoint.x() * AScale, y = APoint.y() * AScale, z = APoint.z() * AScale;
        double t = y * c - z * s;
        z = y * s + z * c;
        y = t; // about x
        t = z * c - x * s;
        x = z * s + x * c;
        z = t; // about y
        t = x * c - y * s;
        y = x * s + y * c;
        x = t; // about z
        return Point3d(x, y, z) + AShift;
    }
} // namespace

using Blk = Blocking<FacetedGeometry>;

TEST_CASE("hex_scaled_jacobian_is_1_for_a_cube", "[BlockTestSuite]") {
    REQUIRE(Blk::hex_scaled_jacobian(unit_cube()) == Approx(1.0));
}

TEST_CASE("hex_scaled_jacobian_ignores_size_position_and_orientation", "[BlockTestSuite]") {
    std::array<Point3d, 8> moved{};
    const auto cube = unit_cube();
    for (std::size_t c = 0; c < 8; ++c) {
        moved[c] = transformed(cube[c], 7.5, 0.7, Vector3d(-3.0, 12.0, 0.25));
    }
    REQUIRE(Blk::hex_scaled_jacobian(moved) == Approx(1.0));
}

TEST_CASE("hex_scaled_jacobian_tends_to_0_as_a_block_shears_flat", "[BlockTestSuite]") {
    // The top face slid a whole edge-length sideways and brought down to almost nothing: the 3
    // edges at a corner end up very nearly coplanar, which is the distortion this measures.
    auto sheared = unit_cube();
    for (std::size_t c = 4; c < 8; ++c) {
        sheared[c] = Point3d(sheared[c].x() + 1.0, sheared[c].y(), 0.01);
    }
    const double q = Blk::hex_scaled_jacobian(sheared);
    REQUIRE(q > 0.0);
    REQUIRE(q < 0.02);
}

TEST_CASE("hex_scaled_jacobian_is_1_for_a_thin_box", "[BlockTestSuite]") {
    // Deliberately recorded, because it is the property people expect this measure not to have: a
    // sliver 1000 times wider than it is tall is *badly proportioned*, not badly *shaped*, and the
    // normalization by each edge's own length is exactly what makes the difference. Aspect ratio is
    // a separate question, and it needs a separate measure to ask it.
    auto thin = unit_cube();
    for (std::size_t c = 4; c < 8; ++c) {
        thin[c] = Point3d(thin[c].x(), thin[c].y(), 0.001);
    }
    REQUIRE(Blk::hex_scaled_jacobian(thin) == Approx(1.0));
}

TEST_CASE("hex_scaled_jacobian_is_negative_for_a_corner_turned_inside_out", "[BlockTestSuite]") {
    auto folded = unit_cube();
    folded[6] = Point3d(0.2, 0.2, 0.2); // pulled through the block, past its 3 neighbours
    REQUIRE(Blk::hex_scaled_jacobian(folded) < 0.0);
}

TEST_CASE("hex_scaled_jacobian_is_0_for_a_collapsed_edge", "[BlockTestSuite]") {
    auto collapsed = unit_cube();
    collapsed[6] = collapsed[7];
    REQUIRE(Blk::hex_scaled_jacobian(collapsed) == Approx(0.0).margin(1e-15));
}

TEST_CASE("quad_scaled_jacobian_is_1_for_a_square_in_any_plane", "[BlockTestSuite]") {
    REQUIRE(Blk::quad_scaled_jacobian(unit_square()) == Approx(1.0));

    std::array<Point3d, 4> oblique{};
    const auto square = unit_square();
    for (std::size_t c = 0; c < 4; ++c) {
        oblique[c] = transformed(square[c], 4.0, 0.9, Vector3d(1.0, -2.0, 5.0));
    }
    REQUIRE(Blk::quad_scaled_jacobian(oblique) == Approx(1.0));
}

TEST_CASE("quad_scaled_jacobian_reads_the_same_from_a_mirrored_frame", "[BlockTestSuite]") {
    // A face's frame is only settled up to a reflection (see frame_of(Face)), so the same quad can
    // legitimately arrive with its corners walked the other way round. Measured against the quad's
    // own normal rather than a fixed axis, that has to make no difference.
    const auto square = unit_square();
    const std::array<Point3d, 4> mirrored{square[0], square[3], square[2], square[1]};
    REQUIRE(Blk::quad_scaled_jacobian(mirrored) == Approx(Blk::quad_scaled_jacobian(square)));
}

TEST_CASE("quad_scaled_jacobian_is_negative_for_a_reflex_corner", "[BlockTestSuite]") {
    std::array<Point3d, 4> reflex = unit_square();
    reflex[2] = Point3d(0.2, 0.2, 0.0); // pulled back inside, leaving the quad non-convex there
    REQUIRE(Blk::quad_scaled_jacobian(reflex) < 0.0);
}

TEST_CASE("quad_scaled_jacobian_is_0_for_a_collapsed_edge", "[BlockTestSuite]") {
    std::array<Point3d, 4> collapsed = unit_square();
    collapsed[2] = collapsed[1];
    REQUIRE(Blk::quad_scaled_jacobian(collapsed) == Approx(0.0).margin(1e-15));
}

TEST_CASE("hex_mean_ratio_is_1_for_a_cube_and_for_nothing_else", "[BlockTestSuite]") {
    REQUIRE(Blk::hex_mean_ratio(unit_cube()) == Approx(1.0));

    std::array<Point3d, 8> moved{};
    const auto cube = unit_cube();
    for (std::size_t c = 0; c < 8; ++c) {
        moved[c] = transformed(cube[c], 7.5, 0.7, Vector3d(-3.0, 12.0, 0.25));
    }
    REQUIRE(Blk::hex_mean_ratio(moved) == Approx(1.0));
}

TEST_CASE("hex_mean_ratio_sees_the_thin_box_the_scaled_jacobian_calls_perfect", "[BlockTestSuite]") {
    // The whole reason both measures exist. The scaled Jacobian scores this sliver 1 — its corners
    // are all square — and a smoother maximizing it would happily produce a blocking of them.
    auto thin = unit_cube();
    for (std::size_t c = 4; c < 8; ++c) {
        thin[c] = Point3d(thin[c].x(), thin[c].y(), 0.001);
    }
    REQUIRE(Blk::hex_scaled_jacobian(thin) == Approx(1.0));
    REQUIRE(Blk::hex_mean_ratio(thin) < 0.02);
    REQUIRE(Blk::hex_mean_ratio(thin) > 0.0);
}

TEST_CASE("hex_mean_ratio_is_negative_for_a_corner_turned_inside_out", "[BlockTestSuite]") {
    auto folded = unit_cube();
    folded[6] = Point3d(0.2, 0.2, 0.2);
    REQUIRE(Blk::hex_mean_ratio(folded) < 0.0);
}

TEST_CASE("hex_mean_ratio_is_0_for_a_collapsed_edge", "[BlockTestSuite]") {
    auto collapsed = unit_cube();
    collapsed[6] = collapsed[7];
    REQUIRE(Blk::hex_mean_ratio(collapsed) == Approx(0.0).margin(1e-15));
}

TEST_CASE("quad_mean_ratio_is_1_for_a_square_and_falls_off_with_the_aspect_ratio", "[BlockTestSuite]") {
    REQUIRE(Blk::quad_mean_ratio(unit_square()) == Approx(1.0));

    std::array<Point3d, 4> oblique{};
    const auto square = unit_square();
    for (std::size_t c = 0; c < 4; ++c) {
        oblique[c] = transformed(square[c], 4.0, 0.9, Vector3d(1.0, -2.0, 5.0));
    }
    REQUIRE(Blk::quad_mean_ratio(oblique) == Approx(1.0));

    const std::array<Point3d, 4> sliver{
        Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 0.01, 0.0), Point3d(0.0, 0.01, 0.0)};
    REQUIRE(Blk::quad_scaled_jacobian(sliver) == Approx(1.0));
    REQUIRE(Blk::quad_mean_ratio(sliver) == Approx(0.02).epsilon(0.01));
}

TEST_CASE("quad_mean_ratio_is_negative_for_a_reflex_corner", "[BlockTestSuite]") {
    std::array<Point3d, 4> reflex = unit_square();
    reflex[2] = Point3d(0.2, 0.2, 0.0);
    REQUIRE(Blk::quad_mean_ratio(reflex) < 0.0);
}

TEST_CASE("block_quality_reports_the_mean_ratio_in_the_blocks_own_frame", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blk blocking(geom);
    const auto block = blocking.create_hex_block(unit_cube());

    REQUIRE(blocking.block_quality(block) == Approx(1.0));
}

TEST_CASE("block_quality_follows_a_corner_that_moved", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blk blocking(geom);
    const auto block = blocking.create_hex_block(unit_cube());

    auto &map = blocking.cmap();
    for (auto it = map.attributes<0>().begin(), itend = map.attributes<0>().end(); it != itend; ++it) {
        if (it->info().point == Point3d(1.0, 1.0, 1.0)) it->info().point = Point3d(0.2, 0.2, 0.2);
    }

    // Read off the corner positions, so a corner written straight into the map — no move_node(), no
    // geometry rebuilt — is already accounted for. What a smoothing iteration needs.
    REQUIRE(blocking.block_quality(block) < 0.0);
}

TEST_CASE("face_quality_reports_the_mean_ratio_in_the_faces_own_frame", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blk blocking(geom);
    const auto face = blocking.create_quad_block(unit_square());

    REQUIRE(blocking.face_quality(face) == Approx(1.0));
}
