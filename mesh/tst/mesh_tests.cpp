#include <gecko/mesh/UnstructuredMesh.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

using Catch::Approx;
using namespace gecko;

TEST_CASE("SimplicialMesh_Nodes_Variables", "[Unstructured Mesh]") {
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
TEST_CASE("SimplicialMesh_Face_Node_Connectivity", "[Unstructured Mesh]") {
    SimplicialMesh mesh;
    REQUIRE(mesh.nb_nodes() == 0);
    auto n0 = mesh.add_node(0, 0, 0);
    ;
    REQUIRE(mesh.nb_nodes() == 1);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(0, 1, 0);
    REQUIRE(mesh.nb_nodes() == 3);
    auto f1 = mesh.add_face(n0, n1, n2);

    REQUIRE(mesh.nb_faces() == 1);

    auto nodes = mesh.face_nodes(f1);
    REQUIRE(nodes.size() == 3);
    REQUIRE(nodes[0] == n0);
    REQUIRE(nodes[1] == n1);
    REQUIRE(nodes[2] == n2);
}

TEST_CASE("SimplicialMesh_Edge_Node_Connectivity", "[Unstructured Mesh]") {
    SimplicialMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);

    auto e0 = mesh.add_edge(n0, n1);
    REQUIRE(mesh.nb_edges() == 1);

    auto nodes = mesh.edge_nodes(e0);
    REQUIRE(nodes.size() == 2);
    REQUIRE(nodes[0] == n0);
    REQUIRE(nodes[1] == n1);
}

TEST_CASE("SimplicialMesh_Cell_Node_Connectivity", "[Unstructured Mesh]") {
    SimplicialMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(0, 1, 0);
    auto n3 = mesh.add_node(0, 0, 1);

    auto c0 = mesh.add_cell(n0, n1, n2, n3);
    REQUIRE(mesh.nb_cells() == 1);

    auto nodes = mesh.cell_nodes(c0);
    REQUIRE(nodes.size() == 4);
    REQUIRE(nodes[0] == n0);
    REQUIRE(nodes[1] == n1);
    REQUIRE(nodes[2] == n2);
    REQUIRE(nodes[3] == n3);
}

TEST_CASE("NeighborEncoding", "[Unstructured Mesh]") {
    SECTION("Same-dimension neighbor (non-negative index)") {
        Int code_3d = 0; // index of the first element of same dimension
        auto ref_3d = NeighborRef::decode(code_3d);

        REQUIRE(ref_3d.type == NeighborRef::Type::Cell3D);
        REQUIRE(ref_3d.index == 0);

        Int code_3d_large = 42;
        auto ref_3d_large = NeighborRef::decode(code_3d_large);

        REQUIRE(ref_3d_large.type == NeighborRef::Type::Cell3D);
        REQUIRE(ref_3d_large.index == 42);
    }

    SECTION("Lower-dimension boundary entity (negative index)") {
        Int idx = 0; // Premier élément de bord
        Int code = NeighborRef::encode_boundary(idx);

        REQUIRE(code == -1); // Indice signifié pour 0 : -0 - 1 = -1

        auto ref = NeighborRef::decode(code);
        REQUIRE(ref.type == NeighborRef::Type::Boundary2D);
        REQUIRE(ref.index == 0);
    }

    SECTION("No neighbor and no materialized lower-dimension entity") {
        auto ref = NeighborRef::decode(NeighborRef::NO_NEIGHBOR);
        REQUIRE(ref.type == NeighborRef::Type::None);
    }
}

TEST_CASE("SimplicialMesh_BuildConnectivity_IsolatedCell", "[Unstructured Mesh]") {
    SimplicialMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(0, 1, 0);
    auto n3 = mesh.add_node(0, 0, 1);
    mesh.add_cell(n0, n1, n2, n3);

    mesh.build_connectivity();

    auto c0 = CellId{0};
    for (SizeT lf = 0; lf < SimplicialMesh::Cell::Faces; ++lf) {
        REQUIRE(mesh.cell_neighbor(c0, lf).type == NeighborRef::Type::None);
    }
}

TEST_CASE("SimplicialMesh_BuildConnectivity_TwoTets", "[Unstructured Mesh]") {
    SimplicialMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(0, 1, 0);
    auto n3 = mesh.add_node(0, 0, 1);
    auto n4 = mesh.add_node(1, 1, 1);

    // Both tets share the face {n0, n1, n2} (local face 3 for both).
    auto cA = mesh.add_cell(n0, n1, n2, n3);
    auto cB = mesh.add_cell(n0, n1, n2, n4);

    // One of tet A's other (genuinely boundary) faces is explicitly materialized: local face 1 = {n0, n2, n3}.
    auto boundaryFace = mesh.add_face(n0, n2, n3);

    mesh.build_connectivity();

    // Shared internal face: same-dimension neighbor on both sides.
    auto sharedA = mesh.cell_neighbor(cA, 3);
    REQUIRE(sharedA.type == NeighborRef::Type::Cell3D);
    REQUIRE(sharedA.index == static_cast<Int>(cB.value));

    auto sharedB = mesh.cell_neighbor(cB, 3);
    REQUIRE(sharedB.type == NeighborRef::Type::Cell3D);
    REQUIRE(sharedB.index == static_cast<Int>(cA.value));

    // Boundary face with a materialized Face: lower-dimension neighbor.
    auto boundary = mesh.cell_neighbor(cA, 1);
    REQUIRE(boundary.type == NeighborRef::Type::Boundary2D);
    REQUIRE(boundary.index == static_cast<Int>(boundaryFace.value));

    auto cells = mesh.face_cells(boundaryFace);
    REQUIRE(cells[0] == cA);
    REQUIRE(cells[1] == CellId::invalid_id);

    // Other boundary faces of tet A were never materialized: no neighbor at all.
    REQUIRE(mesh.cell_neighbor(cA, 0).type == NeighborRef::Type::None);
    REQUIRE(mesh.cell_neighbor(cA, 2).type == NeighborRef::Type::None);
}

TEST_CASE("SimplicialMesh_BuildConnectivity_NonManifold_Throws", "[Unstructured Mesh]") {
    SimplicialMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(0, 1, 0);
    auto n3 = mesh.add_node(0, 0, 1);
    auto n4 = mesh.add_node(1, 1, 1);
    auto n5 = mesh.add_node(-1, -1, -1);

    // Three tets all sharing the same facet {n0, n1, n2} (local face 3 for each): non-manifold.
    mesh.add_cell(n0, n1, n2, n3);
    mesh.add_cell(n0, n1, n2, n4);
    mesh.add_cell(n0, n1, n2, n5);

    REQUIRE_THROWS_AS(mesh.build_connectivity(), std::logic_error);
}

TEST_CASE("SimplicialMesh_BuildConnectivity_FaceToFace", "[Unstructured Mesh]") {
    SimplicialMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(1, 1, 0);
    auto n3 = mesh.add_node(0, 1, 0);

    // Two triangles sharing the edge {n1, n2}.
    auto fA = mesh.add_face(n0, n1, n2); // local edges: {n0,n1}, {n1,n2}, {n2,n0}
    auto fB = mesh.add_face(n1, n2, n3); // local edges: {n1,n2}, {n2,n3}, {n3,n1}

    // One of A's other (genuinely free) edges is explicitly materialized.
    auto boundaryEdge = mesh.add_edge(n0, n1);

    mesh.build_connectivity();

    auto shared = mesh.face_neighbor(fA, 1);
    REQUIRE(shared.type == NeighborRef::Type::Cell3D);
    REQUIRE(shared.index == static_cast<Int>(fB.value));

    auto sharedBack = mesh.face_neighbor(fB, 0);
    REQUIRE(sharedBack.type == NeighborRef::Type::Cell3D);
    REQUIRE(sharedBack.index == static_cast<Int>(fA.value));

    auto boundary = mesh.face_neighbor(fA, 0);
    REQUIRE(boundary.type == NeighborRef::Type::Boundary2D);
    REQUIRE(boundary.index == static_cast<Int>(boundaryEdge.value));

    REQUIRE(mesh.face_neighbor(fA, 2).type == NeighborRef::Type::None);
    REQUIRE(mesh.face_neighbor(fB, 1).type == NeighborRef::Type::None);
    REQUIRE(mesh.face_neighbor(fB, 2).type == NeighborRef::Type::None);
}

TEST_CASE("CubicMesh_BuildConnectivity_TwoHexes", "[Unstructured Mesh]") {
    CubicMesh mesh;
    // Hex A: unit cube, bottom (a0,a1,a2,a3) and top (a4,a5,a6,a7) directly above.
    auto a0 = mesh.add_node(0, 0, 0);
    auto a1 = mesh.add_node(1, 0, 0);
    auto a2 = mesh.add_node(1, 1, 0);
    auto a3 = mesh.add_node(0, 1, 0);
    auto a4 = mesh.add_node(0, 0, 1);
    auto a5 = mesh.add_node(1, 0, 1);
    auto a6 = mesh.add_node(1, 1, 1);
    auto a7 = mesh.add_node(0, 1, 1);
    auto cA = mesh.add_cell(a0, a1, a2, a3, a4, a5, a6, a7);

    // Hex B extends further along +x, reusing A's "right" face {a1,a2,a6,a5} as its "left" face
    // (local {3,0,4,7}): b3=a1, b0=a2, b4=a6, b7=a5.
    auto b1 = mesh.add_node(2, 0, 0);
    auto b2 = mesh.add_node(2, 1, 0);
    auto b5 = mesh.add_node(2, 0, 1);
    auto b6 = mesh.add_node(2, 1, 1);
    auto cB = mesh.add_cell(a2, b1, b2, a1, a6, b5, b6, a5);

    mesh.build_connectivity();

    auto right = mesh.cell_neighbor(cA, 3); // right
    REQUIRE(right.type == NeighborRef::Type::Cell3D);
    REQUIRE(right.index == static_cast<Int>(cB.value));

    auto left = mesh.cell_neighbor(cB, 5); // left
    REQUIRE(left.type == NeighborRef::Type::Cell3D);
    REQUIRE(left.index == static_cast<Int>(cA.value));

    // Every other local face of both hexes is a genuine, un-materialized boundary.
    for (SizeT lf = 0; lf < CubicTraits::Cell::Faces; ++lf) {
        if (lf != 3) {
            REQUIRE(mesh.cell_neighbor(cA, lf).type == NeighborRef::Type::None);
        }
        if (lf != 5) {
            REQUIRE(mesh.cell_neighbor(cB, lf).type == NeighborRef::Type::None);
        }
    }
}
