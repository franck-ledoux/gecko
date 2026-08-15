#include <filesystem>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <gecko/io/GmshCommon.h>
#include <gecko/io/GmshMeshReader.h>
#include <gecko/io/GmshMeshWriter.h>

using Catch::Approx;
using namespace gecko;

namespace {
    /** @brief Path used to write/read back a round-trip test's mesh. */
    std::filesystem::path roundtrip_path(const std::string &name) {
        return std::filesystem::temp_directory_path() / name;
    }
} // namespace

TEST_CASE("Gmsh_Simplicial_WriteRead_Roundtrip", "[Gmsh IO]") {
    SimplicialMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(0, 1, 0);
    auto n3 = mesh.add_node(0, 0, 1);

    GroupRegistry groups;
    auto ridge = groups.add_group("Ridge", GroupDim::Dim1);
    auto wall = groups.add_group("Wall", GroupDim::Dim2);
    auto fluid = groups.add_group("Fluid", GroupDim::Dim3);

    auto &edge_group = mesh.add_variable<GroupId, CellType::Edge>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    auto &edge_entity = mesh.add_variable<Int, CellType::Edge>(std::string(io::ENTITY_TAG_VARIABLE));
    auto &face_group = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    auto &face_entity = mesh.add_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
    auto &cell_group = mesh.add_variable<GroupId, CellType::Cell>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    // Cell entity tag is intentionally left unset, to exercise the "physical only" fallback.

    auto e0 = mesh.add_edge(n0, n1);
    edge_group[e0.value] = ridge;
    edge_entity[e0.value] = 100;

    auto f0 = mesh.add_face(n0, n1, n2);
    face_group[f0.value] = wall;
    face_entity[f0.value] = 7;

    auto c0 = mesh.add_cell(n0, n1, n2, n3);
    cell_group[c0.value] = fluid;

    const auto path = roundtrip_path("gecko_gmsh_simplicial_roundtrip_test.msh").string();
    io::SimplicialMeshWriter::write(path, mesh, groups);

    auto data = io::SimplicialMeshReader::read(path);
    std::filesystem::remove(path);

    REQUIRE(data.mesh.nb_nodes() == 4);
    REQUIRE(data.mesh.nb_edges() == 1);
    REQUIRE(data.mesh.nb_faces() == 1);
    REQUIRE(data.mesh.nb_cells() == 1);

    REQUIRE(data.mesh.node(NodeId{0}).x() == Approx(0.0));
    REQUIRE(data.mesh.node(NodeId{1}).x() == Approx(1.0));
    REQUIRE(data.mesh.node(NodeId{2}).y() == Approx(1.0));
    REQUIRE(data.mesh.node(NodeId{3}).z() == Approx(1.0));

    auto edge_nodes = data.mesh.edge_nodes(EdgeId{0});
    REQUIRE(edge_nodes[0] == NodeId{0});
    REQUIRE(edge_nodes[1] == NodeId{1});

    auto face_nodes = data.mesh.face_nodes(FaceId{0});
    REQUIRE(face_nodes[0] == NodeId{0});
    REQUIRE(face_nodes[1] == NodeId{1});
    REQUIRE(face_nodes[2] == NodeId{2});

    auto cell_nodes = data.mesh.cell_nodes(CellId{0});
    REQUIRE(cell_nodes[0] == NodeId{0});
    REQUIRE(cell_nodes[1] == NodeId{1});
    REQUIRE(cell_nodes[2] == NodeId{2});
    REQUIRE(cell_nodes[3] == NodeId{3});

    REQUIRE(data.groups.size() == 3);
    REQUIRE(data.groups.get_group(GroupId{0}).name == "Ridge");
    REQUIRE(data.groups.get_group(GroupId{0}).dimension == GroupDim::Dim1);
    REQUIRE(data.groups.get_group(GroupId{1}).name == "Wall");
    REQUIRE(data.groups.get_group(GroupId{1}).dimension == GroupDim::Dim2);
    REQUIRE(data.groups.get_group(GroupId{2}).name == "Fluid");
    REQUIRE(data.groups.get_group(GroupId{2}).dimension == GroupDim::Dim3);

    const auto &read_edge_group =
        data.mesh.get_variable<GroupId, CellType::Edge>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    REQUIRE(read_edge_group[0] == GroupId{0});
    const auto &read_edge_entity = data.mesh.get_variable<Int, CellType::Edge>(std::string(io::ENTITY_TAG_VARIABLE));
    REQUIRE(read_edge_entity[0] == 100);

    const auto &read_face_group =
        data.mesh.get_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    REQUIRE(read_face_group[0] == GroupId{1});
    const auto &read_face_entity = data.mesh.get_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
    REQUIRE(read_face_entity[0] == 7);

    const auto &read_cell_group =
        data.mesh.get_variable<GroupId, CellType::Cell>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    REQUIRE(read_cell_group[0] == GroupId{2});
    // No entity tag was set for the cell, so the writer fell back to the physical tag for both,
    // and the reader reads that fallback value back as the (only) entity tag it saw.
    const auto &read_cell_entity = data.mesh.get_variable<Int, CellType::Cell>(std::string(io::ENTITY_TAG_VARIABLE));
    REQUIRE(read_cell_entity[0] == static_cast<Int>(GroupId{2}.value + 1));
}

TEST_CASE("Gmsh_Cubic_WriteRead_Roundtrip", "[Gmsh IO]") {
    CubicMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(1, 1, 0);
    auto n3 = mesh.add_node(0, 1, 0);
    auto n4 = mesh.add_node(0, 0, 1);
    auto n5 = mesh.add_node(1, 0, 1);
    auto n6 = mesh.add_node(1, 1, 1);
    auto n7 = mesh.add_node(0, 1, 1);

    GroupRegistry groups;
    auto bottom = groups.add_group("Bottom", GroupDim::Dim2);

    auto &face_group = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));

    auto f0 = mesh.add_face(n0, n1, n2, n3);
    face_group[f0.value] = bottom;

    mesh.add_cell(n0, n1, n2, n3, n4, n5, n6, n7);

    const auto path = roundtrip_path("gecko_gmsh_cubic_roundtrip_test.msh").string();
    io::CubicMeshWriter::write(path, mesh, groups);

    auto data = io::CubicMeshReader::read(path);
    std::filesystem::remove(path);

    REQUIRE(data.mesh.nb_nodes() == 8);
    REQUIRE(data.mesh.nb_faces() == 1);
    REQUIRE(data.mesh.nb_cells() == 1);

    auto face_nodes = data.mesh.face_nodes(FaceId{0});
    REQUIRE(face_nodes.size() == 4);
    REQUIRE(face_nodes[0] == NodeId{0});
    REQUIRE(face_nodes[3] == NodeId{3});

    auto cell_nodes = data.mesh.cell_nodes(CellId{0});
    REQUIRE(cell_nodes.size() == 8);
    REQUIRE(cell_nodes[7] == NodeId{7});

    REQUIRE(data.groups.size() == 1);
    REQUIRE(data.groups.get_group(GroupId{0}).name == "Bottom");

    const auto &read_face_group =
        data.mesh.get_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    REQUIRE(read_face_group[0] == GroupId{0});
}

TEST_CASE("Gmsh_Read_MissingFile_Throws", "[Gmsh IO]") {
    REQUIRE_THROWS_AS(io::SimplicialMeshReader::read("/nonexistent/path/to/a/file.msh"), std::runtime_error);
}
