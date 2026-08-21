#include <array>
#include <filesystem>

#include <fstream>
#include <sstream>
#include <string>

#include <gecko/block/Blocking.h>
#include <gecko/geom/FacetedGeometry.h>
#include <gecko/io/GmshMeshWriter.h>
#include <gecko/io/VtkMeshWriter.h>
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

        const auto path = (std::filesystem::temp_directory_path() / "gecko_block_mesh_gen_test.msh").string();
        io::SimplicialMeshWriter::write(path, mesh, groups);
        FacetedGeometry geom(path);
        std::filesystem::remove(path);
        return geom;
    }
} // namespace

TEST_CASE("single_quad_block_to_mesh_at_s1_reproduces_corners", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    const std::array<Point3d, 4> corners = {
        Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)};
    blocking.create_quad_block(corners);

    auto mesh = blocking.to_mesh(1);
    REQUIRE(mesh.nb_nodes() == 4);
    REQUIRE(mesh.nb_faces() == 1);
    REQUIRE(mesh.nb_cells() == 0);
    for (UInt i = 0; i < mesh.nb_nodes(); ++i) {
        bool matches_a_corner = false;
        for (const auto &c : corners) {
            if (mesh.node(NodeId{i}) == c) matches_a_corner = true;
        }
        REQUIRE(matches_a_corner);
    }
}

TEST_CASE("single_hex_block_to_mesh_at_s1_reproduces_corners", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    const std::array<Point3d, 8> corners = {Point3d(0.0, 0.0, 0.0),
                                            Point3d(1.0, 0.0, 0.0),
                                            Point3d(1.0, 1.0, 0.0),
                                            Point3d(0.0, 1.0, 0.0),
                                            Point3d(0.0, 0.0, 1.0),
                                            Point3d(1.0, 0.0, 1.0),
                                            Point3d(1.0, 1.0, 1.0),
                                            Point3d(0.0, 1.0, 1.0)};
    blocking.create_hex_block(corners);

    auto mesh = blocking.to_mesh(1);
    REQUIRE(mesh.nb_nodes() == 8);
    REQUIRE(mesh.nb_faces() == 0);
    REQUIRE(mesh.nb_cells() == 1);
    for (UInt i = 0; i < mesh.nb_nodes(); ++i) {
        bool matches_a_corner = false;
        for (const auto &c : corners) {
            if (mesh.node(NodeId{i}) == c) matches_a_corner = true;
        }
        REQUIRE(matches_a_corner);
    }
}

TEST_CASE("single_quad_block_to_mesh_subdivided_has_expected_counts_and_geometry", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    const std::array<Point3d, 4> corners = {
        Point3d(0.0, 0.0, 0.0), Point3d(2.0, 0.0, 0.0), Point3d(2.0, 2.0, 0.0), Point3d(0.0, 2.0, 0.0)};
    blocking.create_quad_block(corners);

    const SizeT s = 3;
    auto mesh = blocking.to_mesh(s);
    REQUIRE(mesh.nb_nodes() == (s + 1) * (s + 1));
    REQUIRE(mesh.nb_faces() == s * s);
    REQUIRE(mesh.nb_cells() == 0);

    // Every generated node must lie within the block's bounding box (a straight quad, so this is a
    // simple, representation-agnostic sanity check that no sample point is wildly wrong).
    for (UInt i = 0; i < mesh.nb_nodes(); ++i) {
        const Point3d &p = mesh.node(NodeId{i});
        REQUIRE(p.x() >= Approx(0.0).margin(1e-9));
        REQUIRE(p.x() <= Approx(2.0).margin(1e-9));
        REQUIRE(p.y() >= Approx(0.0).margin(1e-9));
        REQUIRE(p.y() <= Approx(2.0).margin(1e-9));
    }

    // Every quad face's 4 nodes must be planar-consistent (z=0) and form a non-degenerate cell —
    // spot-check via the face-nodes span's size.
    for (UInt f = 0; f < mesh.nb_faces(); ++f) {
        REQUIRE(mesh.face_nodes(FaceId{f}).size() == 4);
    }
}

TEST_CASE("single_hex_block_to_mesh_subdivided_has_expected_counts", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    const std::array<Point3d, 8> corners = {Point3d(0.0, 0.0, 0.0),
                                            Point3d(1.0, 0.0, 0.0),
                                            Point3d(1.0, 1.0, 0.0),
                                            Point3d(0.0, 1.0, 0.0),
                                            Point3d(0.0, 0.0, 1.0),
                                            Point3d(1.0, 0.0, 1.0),
                                            Point3d(1.0, 1.0, 1.0),
                                            Point3d(0.0, 1.0, 1.0)};
    blocking.create_hex_block(corners);

    const SizeT s = 2;
    auto mesh = blocking.to_mesh(s);
    REQUIRE(mesh.nb_nodes() == (s + 1) * (s + 1) * (s + 1));
    REQUIRE(mesh.nb_cells() == s * s * s);
    REQUIRE(mesh.nb_faces() == 0);
    for (UInt c = 0; c < mesh.nb_cells(); ++c) {
        REQUIRE(mesh.cell_nodes(CellId{c}).size() == 8);
    }
}

TEST_CASE("two_sewn_hex_blocks_to_mesh_share_seam_nodes_without_duplication", "[BlockTestSuite]") {
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
    blocking.build_connectivity();
    REQUIRE(blocking.is_valid_topology());

    const SizeT s = 3;
    auto mesh = blocking.to_mesh(s);

    // Each cube alone would have (s+1)^3 nodes; the shared x=1 face's (s+1)^2 nodes must be counted
    // exactly once, not twice — a hard count, not a coordinate-tolerance heuristic.
    const UInt per_cube = static_cast<UInt>((s + 1) * (s + 1) * (s + 1));
    const UInt shared_face = static_cast<UInt>((s + 1) * (s + 1));
    REQUIRE(mesh.nb_nodes() == 2 * per_cube - shared_face);
    REQUIRE(mesh.nb_cells() == 2 * s * s * s);
    REQUIRE(mesh.nb_faces() == 0);
}

TEST_CASE("two_sewn_quad_blocks_to_mesh_share_seam_nodes_without_duplication", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    const std::array<Point3d, 4> corners_a = {
        Point3d(0.0, 0.0, 0.0), Point3d(1.0, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)};
    const std::array<Point3d, 4> corners_b = {
        Point3d(1.0, 0.0, 0.0), Point3d(2.0, 0.0, 0.0), Point3d(2.0, 1.0, 0.0), Point3d(1.0, 1.0, 0.0)};
    blocking.create_quad_block(corners_a);
    blocking.create_quad_block(corners_b);
    blocking.build_connectivity();
    REQUIRE(blocking.is_valid_topology());

    const SizeT s = 4;
    auto mesh = blocking.to_mesh(s);

    const UInt per_quad = static_cast<UInt>((s + 1) * (s + 1));
    const UInt shared_edge = static_cast<UInt>(s + 1);
    REQUIRE(mesh.nb_nodes() == 2 * per_quad - shared_edge);
    REQUIRE(mesh.nb_faces() == 2 * s * s);
    REQUIRE(mesh.nb_cells() == 0);
}

TEST_CASE("mixed_2d_and_3d_blocking_to_mesh_produces_both_faces_and_cells", "[BlockTestSuite]") {
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    const std::array<Point3d, 4> quad_corners = {
        Point3d(10.0, 0.0, 0.0), Point3d(11.0, 0.0, 0.0), Point3d(11.0, 1.0, 0.0), Point3d(10.0, 1.0, 0.0)};
    const std::array<Point3d, 8> hex_corners = {Point3d(0.0, 0.0, 0.0),
                                                Point3d(1.0, 0.0, 0.0),
                                                Point3d(1.0, 1.0, 0.0),
                                                Point3d(0.0, 1.0, 0.0),
                                                Point3d(0.0, 0.0, 1.0),
                                                Point3d(1.0, 0.0, 1.0),
                                                Point3d(1.0, 1.0, 1.0),
                                                Point3d(0.0, 1.0, 1.0)};
    blocking.create_quad_block(quad_corners);
    blocking.create_hex_block(hex_corners);

    const SizeT s = 2;
    auto mesh = blocking.to_mesh(s);
    REQUIRE(mesh.nb_faces() == s * s);
    REQUIRE(mesh.nb_cells() == s * s * s);
}

TEST_CASE("cubic_hex_block_to_mesh_reflects_curved_geometry", "[BlockTestSuite]") {
    // A block whose edges are curved (via classify() bulging the bottom edge, mirroring
    // blocking_classification_tests.cpp) must produce mesh nodes that follow the curve, not a
    // straight blend — i.e. to_mesh() samples the block's *actual* stored geometry, not just its
    // corners.
    SimplicialMesh mesh_src;
    auto v0 = mesh_src.add_node(0, 0, 0);
    auto v1 = mesh_src.add_node(1, 0, 0);
    auto v2 = mesh_src.add_node(1, 1, 0);
    auto v3 = mesh_src.add_node(0, 1, 0);
    auto vm = mesh_src.add_node(0.5, 0.4, 0);
    GroupRegistry groups;
    auto vtx_group = groups.add_group("Vertices", GroupDim::Dim0);
    auto curve_group = groups.add_group("Curves", GroupDim::Dim1);
    auto &node_group = mesh_src.add_variable<GroupId, CellType::Node>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    auto &node_entity = mesh_src.add_variable<Int, CellType::Node>(std::string(io::ENTITY_TAG_VARIABLE));
    node_group[v0.value] = vtx_group;
    node_entity[v0.value] = 1;
    node_group[v1.value] = vtx_group;
    node_entity[v1.value] = 2;
    node_group[v2.value] = vtx_group;
    node_entity[v2.value] = 3;
    node_group[v3.value] = vtx_group;
    node_entity[v3.value] = 4;
    auto &edge_group = mesh_src.add_variable<GroupId, CellType::Edge>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    auto &edge_entity = mesh_src.add_variable<Int, CellType::Edge>(std::string(io::ENTITY_TAG_VARIABLE));
    auto e0 = mesh_src.add_edge(v0, vm);
    edge_group[e0.value] = curve_group;
    edge_entity[e0.value] = 10;
    auto e1 = mesh_src.add_edge(vm, v1);
    edge_group[e1.value] = curve_group;
    edge_entity[e1.value] = 10;
    auto e2 = mesh_src.add_edge(v1, v2);
    edge_group[e2.value] = curve_group;
    edge_entity[e2.value] = 11;
    auto e3 = mesh_src.add_edge(v2, v3);
    edge_group[e3.value] = curve_group;
    edge_entity[e3.value] = 12;
    auto e4 = mesh_src.add_edge(v3, v0);
    edge_group[e4.value] = curve_group;
    edge_entity[e4.value] = 13;

    const auto path = (std::filesystem::temp_directory_path() / "gecko_block_mesh_gen_bent_test.msh").string();
    io::SimplicialMeshWriter::write(path, mesh_src, groups);
    const FacetedGeometry geom(path);
    std::filesystem::remove(path);

    Blocking<FacetedGeometry, BezierCurve<3, Point3d>> blocking(geom);
    const std::array<Point3d, 4> corners = {
        Point3d(0.01, 0.0, 0.0), Point3d(0.99, 0.0, 0.0), Point3d(1.0, 1.0, 0.0), Point3d(0.0, 1.0, 0.0)};
    blocking.create_quad_block(corners);
    blocking.classify(0.05, 0.35);

    const SizeT s = 4;
    auto out_mesh = blocking.to_mesh(s);
    REQUIRE(out_mesh.nb_nodes() == (s + 1) * (s + 1));

    // At least one generated node must lie clearly off the z=0,y=0 straight baseline (proving the
    // bulge reached the actual mesh output), and none may have wildly wrong values.
    bool found_bulged_node = false;
    for (UInt i = 0; i < out_mesh.nb_nodes(); ++i) {
        const Point3d &p = out_mesh.node(NodeId{i});
        if (p.y() > 0.1 && p.x() > 0.0 && p.x() < 1.0) found_bulged_node = true;
    }
    REQUIRE(found_bulged_node);
}

TEST_CASE("quad_block_to_mesh_propagates_node_classification_from_owning_cells", "[BlockTestSuite]") {
    // Classification is assigned by hand rather than through classify()'s geometric snapping: this
    // isolates to_mesh()'s own propagation logic (the thing this test exists for) from classify()'s
    // separate, already-covered snapping behavior (see blocking_classification_tests.cpp).
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    const std::array<Point3d, 4> corners = {
        Point3d(0.0, 0.0, 0.0), Point3d(2.0, 0.0, 0.0), Point3d(2.0, 2.0, 0.0), Point3d(0.0, 2.0, 0.0)};
    auto face = blocking.create_quad_block(corners);

    auto node_it = blocking.cmap().attributes<0>().begin();
    node_it->info().geom_targets = {{GroupDim::Dim0, 7}};

    auto edge_it = blocking.cmap().attributes<1>().begin();
    edge_it->info().geom_targets = {{GroupDim::Dim1, 11}};

    face->info().geom_targets = {{GroupDim::Dim2, 1}};

    const SizeT s = 2;
    auto mesh = blocking.to_mesh(s);
    auto &dims = mesh.get_variable<Int, CellType::Node>(
        std::string(Blocking<FacetedGeometry>::NODE_CLASSIFICATION_DIM_VARIABLE));
    auto &tags = mesh.get_variable<Int, CellType::Node>(
        std::string(Blocking<FacetedGeometry>::NODE_CLASSIFICATION_TAG_VARIABLE));
    REQUIRE(dims.size() == mesh.nb_nodes());
    REQUIRE(tags.size() == mesh.nb_nodes());

    // s=2 gives exactly 1 interior sample point per edge (dim0=1) and 1 for the face (dim0=1,1);
    // the tagged corner/edge/face must each contribute exactly that many matching nodes, and every
    // other node (the 3 untouched corners, 3 untouched edges' interior points) must read -1/-1.
    int classified_as_vertex = 0, classified_as_curve = 0, classified_as_surface = 0, unclassified = 0;
    for (UInt i = 0; i < mesh.nb_nodes(); ++i) {
        const Int dim = dims[i];
        const Int tag = tags[i];
        if (dim == -1) {
            REQUIRE(tag == -1);
            ++unclassified;
        } else if (dim == static_cast<Int>(GroupDim::Dim0)) {
            REQUIRE(tag == 7);
            ++classified_as_vertex;
        } else if (dim == static_cast<Int>(GroupDim::Dim1)) {
            REQUIRE(tag == 11);
            ++classified_as_curve;
        } else {
            REQUIRE(dim == static_cast<Int>(GroupDim::Dim2));
            REQUIRE(tag == 1);
            ++classified_as_surface;
        }
    }
    REQUIRE(classified_as_vertex == 1);
    REQUIRE(classified_as_curve == 1);
    REQUIRE(classified_as_surface == 1);
    REQUIRE(unclassified == static_cast<int>(mesh.nb_nodes()) - 3);
}

TEST_CASE("hex_block_to_mesh_at_s1_can_be_exported_to_vtk_legacy", "[BlockTestSuite]") {
    // Issue #22's "VTK legacy export of the block structure" ask: to_mesh(1) reproduces exactly
    // the coarse block topology (block corners, one hex per top-cell), suitable for VtkMeshWriter.
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    const std::array<Point3d, 8> corners = {Point3d(0.0, 0.0, 0.0),
                                            Point3d(1.0, 0.0, 0.0),
                                            Point3d(1.0, 1.0, 0.0),
                                            Point3d(0.0, 1.0, 0.0),
                                            Point3d(0.0, 0.0, 1.0),
                                            Point3d(1.0, 0.0, 1.0),
                                            Point3d(1.0, 1.0, 1.0),
                                            Point3d(0.0, 1.0, 1.0)};
    blocking.create_hex_block(corners);

    const auto path = (std::filesystem::temp_directory_path() / "gecko_block_to_vtk_test.vtk").string();
    io::CubicVtkWriter::write(path, blocking.to_mesh(1));

    std::ifstream in(path);
    REQUIRE(in.is_open());
    std::ostringstream ss;
    ss << in.rdbuf();
    const std::string content = ss.str();
    // Close before removing: Windows (unlike POSIX) refuses to delete a file with an open handle.
    in.close();
    std::filesystem::remove(path);

    REQUIRE(content.find("DATASET UNSTRUCTURED_GRID") != std::string::npos);
    REQUIRE(content.find("POINTS 8 double") != std::string::npos);
    REQUIRE(content.find("CELLS 1 9") != std::string::npos);
    REQUIRE(content.find("\n12\n") != std::string::npos); // VTK_HEXAHEDRON
}

TEST_CASE("hex_block_to_mesh_classification_variables_can_be_exported_to_vtk_legacy", "[BlockTestSuite]") {
    // The 2 node classification variables to_mesh() always writes are ordinary named node
    // variables from VtkMeshWriter's point of view — this exercises passing their names through to
    // write(), end to end from an unclassified block (everything reads -1/-1).
    const FacetedGeometry geom = make_minimal_geom_model();
    Blocking<FacetedGeometry> blocking(geom);
    const std::array<Point3d, 8> corners = {Point3d(0.0, 0.0, 0.0),
                                            Point3d(1.0, 0.0, 0.0),
                                            Point3d(1.0, 1.0, 0.0),
                                            Point3d(0.0, 1.0, 0.0),
                                            Point3d(0.0, 0.0, 1.0),
                                            Point3d(1.0, 0.0, 1.0),
                                            Point3d(1.0, 1.0, 1.0),
                                            Point3d(0.0, 1.0, 1.0)};
    blocking.create_hex_block(corners);

    const auto path = (std::filesystem::temp_directory_path() / "gecko_block_to_vtk_classification_test.vtk").string();
    io::CubicVtkWriter::write(path,
                              blocking.to_mesh(1),
                              {std::string(Blocking<FacetedGeometry>::NODE_CLASSIFICATION_DIM_VARIABLE),
                               std::string(Blocking<FacetedGeometry>::NODE_CLASSIFICATION_TAG_VARIABLE)});

    std::ifstream in(path);
    REQUIRE(in.is_open());
    std::ostringstream ss;
    ss << in.rdbuf();
    const std::string content = ss.str();
    in.close();
    std::filesystem::remove(path);

    REQUIRE(content.find("POINT_DATA 8") != std::string::npos);
    REQUIRE(content.find("SCALARS classification_dim int 1\nLOOKUP_TABLE default\n-1\n-1\n-1\n-1\n-1\n-1\n-1\n-1\n") !=
            std::string::npos);
    REQUIRE(content.find("SCALARS classification_tag int 1\nLOOKUP_TABLE default\n-1\n-1\n-1\n-1\n-1\n-1\n-1\n-1\n") !=
            std::string::npos);
}
