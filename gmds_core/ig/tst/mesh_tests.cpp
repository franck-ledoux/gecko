#include <catch2/catch_test_macros.hpp>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>

TEST_CASE("Mesh Node operations", "[MeshNode]") {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
    gmds::Node n0 = m.newNode(0, 0, 0);
    REQUIRE(n0.point() == gmds::math::Point(0, 0, 0));

    n0.setX(1);
    REQUIRE(n0.X() == 1);

    n0.setY(2);
    REQUIRE(n0.Y() == 2);

    n0.setZ(3);
    REQUIRE(n0.Z() == 3);

    n0.setXYZ(4, 5, 6);
    REQUIRE(n0.point() == gmds::math::Point(4, 5, 6));

    n0.setPoint(gmds::math::Point(1, 1, 1));
    REQUIRE(n0.point() == gmds::math::Point(1, 1, 1));
}

TEST_CASE("Mesh Node and Face Iteration", "[MeshNodeFaceIter]") {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

    gmds::Node n0 = m.newNode(0, 0, 0);
    gmds::Node n1 = m.newNode(1, 1, 1);
    gmds::Node n2 = m.newNode(0, 1, 1);
    gmds::Node n3 = m.newNode(1, 0, 1);

    m.newTriangle(n0, n1, n3);
    m.newTriangle(n0, n3, n2);
    m.newQuad(n0, n1, n2, n3);

    size_t nodeCount = 0;
    for (auto i : m.nodes()) {
        nodeCount++;
    }
    REQUIRE(nodeCount == 4);

    m.deleteNode(2);

    size_t remainingNodes = 0;
    for (auto i : m.nodes()) {
        remainingNodes++;
    }
    REQUIRE(remainingNodes == 3);
}

TEST_CASE("Node Variables and Groups", "[NodeVarGroup]") {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

    gmds::Node n0 = m.newNode(0, 0, 0);
    gmds::Node n1 = m.newNode(1, 1, 1);
    gmds::Node n2 = m.newNode(0, 1, 1);
    gmds::Node n3 = m.newNode(1, 0, 1);

    auto* var = m.newVariable<int, gmds::GMDS_NODE>("int_var");

    for (auto i : m.nodes()) {
        REQUIRE(var->value(i) == 0);
        var->value(i) = i;
    }

    for (auto i : m.nodes()) {
        REQUIRE((*var)[i] == i);
    }

    auto* var2 = m.getVariable<int, gmds::GMDS_NODE>("int_var");

    for (auto i : m.nodes()) {
        REQUIRE((*var2)[i] == i);
    }

    auto* grp = m.newGroup<gmds::Node>("cloud");
    grp->add(n0);
    grp->add(n1);
    grp->add(n2.id());

    REQUIRE(grp->size() == 3);
    REQUIRE((*grp)[0] == n0.id());
    REQUIRE((*grp)[1] == n1.id());
    REQUIRE((*grp)[2] == n2.id());
}

TEST_CASE("Mesh Accessor Validation", "[MeshAccessor]") {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::F2N | gmds::F2E | gmds::E2N | gmds::E2F | gmds::N2E | gmds::N2F));

    REQUIRE(m.getDim() == 3);

    gmds::Node n0 = m.newNode(0, 0, 1);
    gmds::Node n1 = m.newNode(1, 0, 1);
    gmds::Node n2 = m.newNode(1, 1, 1);
    gmds::Node n3 = m.newNode(0, 1, 1);
    gmds::Node n4 = m.newNode(2, 0, 1);
    gmds::Node n5 = m.newNode(2, 1, 1);

    auto f1 = m.newTriangle(n0, n1, n3);
    auto f2 = m.newTriangle(n1, n3, n2);
    auto f3 = m.newQuad(n1, n2, n5, n4);

    gmds::MeshDoctor doc(&m);
    doc.buildEdgesAndX2E();

    REQUIRE(f1.getIDs<gmds::Edge>().size() == 3);
    REQUIRE(f2.getIDs<gmds::Edge>().size() == 3);
    REQUIRE(f3.getIDs<gmds::Edge>().size() == 4);

    doc.updateUpwardConnectivity();

    REQUIRE(f1.getAllIDs<gmds::Node>().size() == 3);
    REQUIRE(f2.getAllIDs<gmds::Node>().size() == 3);
    REQUIRE(f3.getAllIDs<gmds::Node>().size() == 4);
}
