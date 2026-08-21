#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <gecko/io/VtkMeshWriter.h>

using Catch::Approx;
using namespace gecko;

namespace {
    /** @brief Path used to write a VTK export test's mesh. */
    std::filesystem::path vtk_path(const std::string &name) { return std::filesystem::temp_directory_path() / name; }

    /** @brief Reads the whole content of a small text file into a string. */
    std::string read_file(const std::filesystem::path &path) {
        std::ifstream in(path);
        REQUIRE(in.is_open());
        std::ostringstream ss;
        ss << in.rdbuf();
        return ss.str();
    }
} // namespace

TEST_CASE("Vtk_Cubic_SingleHex_HasExpectedHeaderCountsAndConnectivity", "[Vtk IO]") {
    CubicMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(1, 1, 0);
    auto n3 = mesh.add_node(0, 1, 0);
    auto n4 = mesh.add_node(0, 0, 1);
    auto n5 = mesh.add_node(1, 0, 1);
    auto n6 = mesh.add_node(1, 1, 1);
    auto n7 = mesh.add_node(0, 1, 1);
    mesh.add_cell(n0, n1, n2, n3, n4, n5, n6, n7);

    const auto path = vtk_path("gecko_vtk_single_hex_test.vtk").string();
    io::CubicVtkWriter::write(path, mesh);
    const std::string content = read_file(path);
    std::filesystem::remove(path);

    REQUIRE(content.find("DATASET UNSTRUCTURED_GRID") != std::string::npos);
    REQUIRE(content.find("POINTS 8 double") != std::string::npos);
    REQUIRE(content.find("CELLS 1 9") != std::string::npos);
    REQUIRE(content.find("CELL_TYPES 1") != std::string::npos);
    REQUIRE(content.find("8 0 1 2 3 4 5 6 7") != std::string::npos);
    REQUIRE(content.find("\n12\n") != std::string::npos); // VTK_HEXAHEDRON
}

TEST_CASE("Vtk_Cubic_MixedFacesAndCells_HasExpectedCounts", "[Vtk IO]") {
    CubicMesh mesh;
    // A standalone quad face...
    auto q0 = mesh.add_node(0, 0, 0);
    auto q1 = mesh.add_node(1, 0, 0);
    auto q2 = mesh.add_node(1, 1, 0);
    auto q3 = mesh.add_node(0, 1, 0);
    mesh.add_face(q0, q1, q2, q3);
    // ...and a separate hex cell.
    auto h0 = mesh.add_node(10, 0, 0);
    auto h1 = mesh.add_node(11, 0, 0);
    auto h2 = mesh.add_node(11, 1, 0);
    auto h3 = mesh.add_node(10, 1, 0);
    auto h4 = mesh.add_node(10, 0, 1);
    auto h5 = mesh.add_node(11, 0, 1);
    auto h6 = mesh.add_node(11, 1, 1);
    auto h7 = mesh.add_node(10, 1, 1);
    mesh.add_cell(h0, h1, h2, h3, h4, h5, h6, h7);

    const auto path = vtk_path("gecko_vtk_mixed_test.vtk").string();
    io::CubicVtkWriter::write(path, mesh);
    const std::string content = read_file(path);
    std::filesystem::remove(path);

    REQUIRE(content.find("POINTS 12 double") != std::string::npos);
    // 1 quad (5 ints: 4 + count) + 1 hex (9 ints: 8 + count) = 14 total.
    REQUIRE(content.find("CELLS 2 14") != std::string::npos);
    REQUIRE(content.find("CELL_TYPES 2") != std::string::npos);
    REQUIRE(content.find("\n9\n") != std::string::npos);  // VTK_QUAD (the face)
    REQUIRE(content.find("\n12\n") != std::string::npos); // VTK_HEXAHEDRON (the cell)
}

TEST_CASE("Vtk_Simplicial_UsesTriangleAndTetraTypes", "[Vtk IO]") {
    SimplicialMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(0, 1, 0);
    mesh.add_face(n0, n1, n2);

    const auto path = vtk_path("gecko_vtk_simplicial_test.vtk").string();
    io::SimplicialVtkWriter::write(path, mesh);
    const std::string content = read_file(path);
    std::filesystem::remove(path);

    REQUIRE(content.find("POINTS 3 double") != std::string::npos);
    REQUIRE(content.find("CELLS 1 4") != std::string::npos); // 1 triangle: 3 + count.
    REQUIRE(content.find("3 0 1 2") != std::string::npos);
    REQUIRE(content.find("\n5\n") != std::string::npos); // VTK_TRIANGLE
}

TEST_CASE("Vtk_EmptyMesh_WritesWellFormedEmptyHeaders", "[Vtk IO]") {
    CubicMesh mesh;
    const auto path = vtk_path("gecko_vtk_empty_test.vtk").string();
    io::CubicVtkWriter::write(path, mesh);
    const std::string content = read_file(path);
    std::filesystem::remove(path);

    REQUIRE(content.find("POINTS 0 double") != std::string::npos);
    REQUIRE(content.find("CELLS 0 0") != std::string::npos);
    REQUIRE(content.find("CELL_TYPES 0") != std::string::npos);
}

TEST_CASE("Vtk_NamedNodeIntVariables_WritesPointDataScalarsSections", "[Vtk IO]") {
    CubicMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(1, 1, 0);
    auto n3 = mesh.add_node(0, 1, 0);
    mesh.add_face(n0, n1, n2, n3);

    auto &dims = mesh.add_variable<Int, CellType::Node>("classification_dim");
    dims[n0.value] = 0;
    dims[n1.value] = -1;
    dims[n2.value] = 2;
    dims[n3.value] = -1;
    auto &tags = mesh.add_variable<Int, CellType::Node>("classification_tag");
    tags[n0.value] = 7;
    tags[n1.value] = -1;
    tags[n2.value] = 1;
    tags[n3.value] = -1;

    const auto path = vtk_path("gecko_vtk_point_data_test.vtk").string();
    io::CubicVtkWriter::write(path, mesh, {"classification_dim", "classification_tag"});
    const std::string content = read_file(path);
    std::filesystem::remove(path);

    REQUIRE(content.find("POINT_DATA 4") != std::string::npos);
    REQUIRE(content.find("SCALARS classification_dim int 1\nLOOKUP_TABLE default\n0\n-1\n2\n-1\n") !=
            std::string::npos);
    REQUIRE(content.find("SCALARS classification_tag int 1\nLOOKUP_TABLE default\n7\n-1\n1\n-1\n") !=
            std::string::npos);
    // Field order must follow the requested list, POINT_DATA before either SCALARS block.
    REQUIRE(content.find("POINT_DATA") < content.find("SCALARS classification_dim"));
    REQUIRE(content.find("SCALARS classification_dim") < content.find("SCALARS classification_tag"));
}

TEST_CASE("Vtk_UnknownNodeIntVariableName_IsSilentlySkipped", "[Vtk IO]") {
    CubicMesh mesh;
    mesh.add_node(0, 0, 0);

    const auto path = vtk_path("gecko_vtk_point_data_missing_test.vtk").string();
    io::CubicVtkWriter::write(path, mesh, {"does_not_exist"});
    const std::string content = read_file(path);
    std::filesystem::remove(path);

    REQUIRE(content.find("POINT_DATA") == std::string::npos);
    REQUIRE(content.find("SCALARS") == std::string::npos);
}
