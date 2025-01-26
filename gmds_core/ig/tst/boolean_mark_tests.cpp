#include <catch2/catch_test_macros.hpp>
#include <gmds/ig/Mesh.h>

TEST_CASE("BooleanMarkMesh: testNodeMarks", "[BooleanMarkMesh]") {
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

	gmds::Node n0 = m.newNode(0, 0, 0);
	gmds::Node n1 = m.newNode(1, 1, 1);
	gmds::Node n2 = m.newNode(0, 1, 1);
	gmds::Node n3 = m.newNode(1, 0, 1);

	m.newTriangle(n0, n1, n3);
	m.newTriangle(n0, n3, n2);
	m.newQuad(n0, n1, n2, n3);

	int mark = m.newMark<gmds::Node>();

	REQUIRE(m.getMarks(n0) == gmds::Marks32());
	REQUIRE(m.getMarks(n1) == gmds::Marks32());
	REQUIRE(m.getMarks(n2) == gmds::Marks32());
	REQUIRE(m.getMarks(n3) == gmds::Marks32());

	int nb_marked_nodes = 0;
	for (auto i : m.nodes()) {
		if (m.isMarked<gmds::Node>(i, mark))
			nb_marked_nodes++;
	}
	REQUIRE(nb_marked_nodes == 0);

	m.mark(n0, mark);
	m.mark<gmds::Node>(n1.id(), mark);
	m.mark<gmds::Node>(n2.id(), mark);

	nb_marked_nodes = 0;
	for (auto i : m.nodes()) {
		if (m.isMarked(m.get<gmds::Node>(i), mark))
			nb_marked_nodes++;
	}
	REQUIRE(nb_marked_nodes == 3);

	gmds::Marks32 marks32;
	marks32.set(mark);
	REQUIRE(m.getMarks(n0) == marks32);
	REQUIRE(m.getMarks(n3) == gmds::Marks32());

	m.negateMaskMark<gmds::Node>(mark);

	nb_marked_nodes = 0;
	for (auto i : m.nodes()) {
		if (m.isMarked(m.get<gmds::Node>(i), mark))
			nb_marked_nodes++;
	}
	REQUIRE(nb_marked_nodes == 1);

	REQUIRE(m.getMarks(n0) == marks32);
	REQUIRE(m.getMarks(n3) == gmds::Marks32());

	m.unmarkAll<gmds::Node>(mark);

	nb_marked_nodes = 0;
	for (auto i : m.nodes()) {
		if (m.isMarked(m.get<gmds::Node>(i), mark))
			nb_marked_nodes++;
	}
	REQUIRE(nb_marked_nodes == 0);

	m.freeMark<gmds::Node>(mark);
}
