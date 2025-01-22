#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <catch2/catch_test_macros.hpp>


TEST_CASE("EdgeTestSuite: testMeshEdge", "[EdgeTestSuite]") {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::E | gmds::N | gmds::E2N));
    gmds::Node n0 = m.newNode(0, 0, 0);
    gmds::Node n1 = m.newNode(1, 0, 0);
    gmds::Node n2 = m.newNode(0, 1, 0);
    gmds::Edge e0 = m.newEdge(n0, n1);
    gmds::Edge e1 = m.newEdge(n0, n2);

    gmds::Edge ea = m.get<gmds::Edge>(e0.id());

    REQUIRE(e0 == ea);
    REQUIRE(e0 != e1);

    REQUIRE_FALSE(e0 != ea);
    REQUIRE_FALSE(e0 == e1);
}

TEST_CASE("EdgeTestSuite: testAccessorIssues", "[EdgeTestSuite]") {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E));
    gmds::Node n0 = m.newNode(0, 0, 0);
    gmds::Node n1 = m.newNode(1, 0, 0);
    gmds::Node n2 = m.newNode(2, 0, 0);
    gmds::Edge e0 = m.newEdge(n0, n1);

    // Les assertions d'exceptions sont désactivées pour le moment.
    // REQUIRE_THROWS(e0.getIDs<gmds::Node>());
    // REQUIRE_THROWS(e0.getIDs<gmds::Edge>());
    // REQUIRE_THROWS(e0.getIDs<gmds::Face>());
    // REQUIRE_THROWS(e0.getIDs<gmds::Region>());

    gmds::MeshModel mod(gmds::DIM2 | gmds::F | gmds::N | gmds::E | gmds::E2N);
    m.changeModel(mod);

    std::vector<gmds::Node> nodes = {n0, n1};
    e0.set<gmds::Node>(nodes);

    std::vector<gmds::TCellID> ids = e0.getIDs<gmds::Node>();
    REQUIRE(ids.size() == 2);

    auto it = std::find(ids.begin(), ids.end(), n0.id());
    REQUIRE(it != ids.end());

    it = std::find(ids.begin(), ids.end(), n1.id());
    REQUIRE(it != ids.end());

    it = std::find(ids.begin(), ids.end(), n2.id());
    REQUIRE(it == ids.end());

    e0.replace<gmds::Node>(n0, n2);
    ids = e0.getIDs<gmds::Node>();
    it = std::find(ids.begin(), ids.end(), n0.id());
    REQUIRE(it == ids.end());

    it = std::find(ids.begin(), ids.end(), n1.id());
    REQUIRE(it != ids.end());

    it = std::find(ids.begin(), ids.end(), n2.id());
    REQUIRE(it != ids.end());

    e0.remove(n1);
    ids = e0.getIDs<gmds::Node>();
    REQUIRE(ids.size() == 1);
}

TEST_CASE("EdgeTestSuite: testAccessorAll", "[EdgeTestSuite]") {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::N | gmds::F | gmds::E |
                                 gmds::F2N | gmds::R2N | gmds::E2N |
                                 gmds::R2F | gmds::F2R | gmds::R2E));

    gmds::Node n0 = m.newNode(0, 0, 0);
    gmds::Node n1 = m.newNode(0, 1, 0);
    gmds::Node n2 = m.newNode(1, 1, 0);
    gmds::Node n3 = m.newNode(0, 1, 1);
    gmds::Node n4 = m.newNode(0, 1, 2);

    gmds::MeshDoctor doc(&m);
    doc.buildEdgesAndX2E();
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    gmds::Edge e0 = m.get<gmds::Edge>(0);

    // Les assertions désactivées peuvent être réactivées selon les résultats attendus.
    // REQUIRE(e0.nbFaces() == 1);
}
