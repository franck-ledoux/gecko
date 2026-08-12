#include <gecko/mesh/UnstructuredMesh.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

using Catch::Approx;
using namespace gecko;

TEST_CASE("UM_Variables", "[Unstructured Mesh]") {
    SimplicialMesh mesh;
    auto &v = mesh.add_variable<int, CellType::Node>("test");
    REQUIRE(v.size() == 0);
    REQUIRE(mesh.has_variable<CellType::Node>("test"));
    mesh.add_node(0, 0, 0);
    REQUIRE(mesh.nb_nodes() == 1);
    REQUIRE(v.size() == 1);
    mesh.remove_variable<CellType::Node>("test");
    REQUIRE_FALSE(mesh.has_variable<CellType::Node>("test"));
}
TEST_CASE("UM_Connectivity", "[Unstructured Mesh]") {
    SimplicialMesh mesh;
    auto &v = mesh.add_variable<int, CellType::Node>("test");
    REQUIRE(mesh.nb_nodes() == 0);
    mesh.add_node(0, 0, 0);
    REQUIRE(v.size() == 1);
    mesh.remove_variable<CellType::Node>("test");
    REQUIRE_FALSE(mesh.has_variable<CellType::Node>("test"));
}
/*
TEST_CASE("Variables", "[Mesh]") {
    SECTION("Tet  (Indices positifs)") {
        int32_t code_3d = 0; // Index du premier élément 3D
        auto ref_3d = decode_neighbor(code_3d);

        REQUIRE(ref_3d.type == NeighborRef::Type::Cell3D);
        REQUIRE(ref_3d.index == 0);

        int32_t code_3d_large = 42;
        auto ref_3d_large = decode_neighbor(code_3d_large);

        REQUIRE(ref_3d_large.type == NeighborRef::Type::Cell3D);
        REQUIRE(ref_3d_large.index == 42);
    }

    SECTION("Triangle/Quad de bord (Indices négatifs)") {
        int32_t tri_idx = 0; // Premier triangle de bord
        int32_t code_2d = encode_boundary(tri_idx);

        REQUIRE(code_2d == -1); // Indice signifié pour 0 : -0 - 1 = -1

        auto ref_2d = decode_neighbor(code_2d);
        REQUIRE(ref_2d.type == NeighborRef::Type::Boundary2D);
        REQUIRE(ref_2d.index == 0);
    }
}*/