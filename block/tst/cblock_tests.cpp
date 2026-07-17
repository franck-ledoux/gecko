#include <gecko/cblock/Blocking.h>
#include <gecko/cblock/BlockingClassifier.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <unit_test_config.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
using Catch::Approx;

void
setUp(gmds::cad::FACManager &AGeomManager)
{
	gmds::Mesh m_vol(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::E | gmds::N |
	                                 gmds::R2N | gmds::R2F | gmds::R2E | gmds::F2N |
	                                 gmds::F2R | gmds::F2E
	                                 | gmds::E2F | gmds::E2N | gmds::N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/tet_in_box.vtk";
	gmds::IGMeshIOService ioService(&m_vol);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::R);
	vtkReader.read(vtk_file);
	gmds::MeshDoctor doc(&m_vol);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	AGeomManager.initFrom3DMesh(&m_vol);
}

std::tuple<int,int,int,int> get_node_statistics(gecko::cblock::Blocking& ABlocking){
    auto nb_on_vertex=0;
    auto nb_on_curve=0;
    auto nb_on_surface=0;
    auto nb_in_volume=0;
    std::vector<gecko::cblock::Blocking::Node> all_nodes = ABlocking.get_all_nodes();
    for(auto n:all_nodes){
        if(n->info().geom_dim==0)
            nb_on_vertex++;
        else if(n->info().geom_dim==1)
            nb_on_curve++;
        else if(n->info().geom_dim==2)
            nb_on_surface++;
        else if(n->info().geom_dim==3)
            nb_in_volume++;
    }
    return std::make_tuple(nb_on_vertex,nb_on_curve,nb_on_surface,nb_in_volume);
}

std::tuple<int,int,int> get_edge_statistics(gecko::cblock::Blocking& ABlocking){
    auto nb_on_curve=0;
    auto nb_on_surface=0;
    auto nb_in_volume=0;
    std::vector<gecko::cblock::Blocking::Edge> all_edges = ABlocking.get_all_edges();
    for(auto e:all_edges){
        if(e->info().geom_dim==1)
            nb_on_curve++;
        else if(e->info().geom_dim==2)
            nb_on_surface++;
        else if(e->info().geom_dim==3)
            nb_in_volume++;
    }
    return std::make_tuple(nb_on_curve,nb_on_surface,nb_in_volume);
}

std::tuple<int,int> get_face_statistics(gecko::cblock::Blocking& ABlocking){
    auto nb_on_surface=0;
    auto nb_in_volume=0;
    std::vector<gecko::cblock::Blocking::Face> all_faces = ABlocking.get_all_faces();
    for(auto f:all_faces){
        if(f->info().geom_dim==2)
            nb_on_surface++;
        else if(f->info().geom_dim==3)
            nb_in_volume++;
    }
    return std::make_tuple(nb_on_surface,nb_in_volume);
}

void export_vtk(gecko::cblock::Blocking& ABlocking, int AModel, const std::string& AFileName){
    gmds::Mesh m_out(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F | gmds::R | gmds::E2N | gmds::F2N | gmds::R2N));
    ABlocking.convert_to_mesh(m_out);
    gmds::IGMeshIOService ioService(&m_out);
    gmds::VTKWriter writer(&ioService);
    writer.setCellOptions(AModel);
    writer.setDataOptions(AModel);
    writer.write(AFileName);
}

TEST_CASE("single_block", "[BlockingTestSuite]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::cblock::Blocking bl(&geom_model);

    gmds::math::Point p000(0, 0, 0);
    gmds::math::Point p010(0, 1, 0);
    gmds::math::Point p110(1, 1, 0);
    gmds::math::Point p100(1, 0, 0);

    gmds::math::Point p001(0, 0, 1);
    gmds::math::Point p011(0, 1, 1);
    gmds::math::Point p111(1, 1, 1);
    gmds::math::Point p101(1, 0, 1);

    auto b = bl.create_block(p000, p010, p110, p100, p001, p011, p111, p101);

    REQUIRE(b->info().geom_dim == 4);
    REQUIRE(b->info().geom_id == gmds::NullID);

    auto fs = bl.get_faces_of_block(b);
    REQUIRE(fs.size() == 6);

    auto block_center = bl.get_center_of_block(b);
    for (auto f : fs) {
        gmds::math::Point face_center = bl.get_center_of_face(f);
        REQUIRE(std::abs(block_center.distance(face_center) - 0.5) < 1e-8);
    }

}

TEST_CASE("BlockingTestSuite - split_one_block_once", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::cblock::Blocking bl(&geom_model, true);
    gecko::cblock::BlockingClassifier cl(&bl);
    cl.classify();

    auto e = bl.get_all_edges()[0];
    auto e_id = e->info().geom_id;
    auto e_dim = e->info().geom_dim;

    bl.cut_sheet(e);
    REQUIRE(bl.get_nb_cells<0>() == 12);
    REQUIRE(bl.get_nb_cells<1>() == 20);
    REQUIRE(bl.get_nb_cells<2>() == 11);
    REQUIRE(bl.get_nb_cells<3>() == 2);
    REQUIRE(bl.cmap()->is_valid());

    int classified_nodes = 0;
    int classified_edges = 0;
    for (auto cur_edge : bl.get_all_edges()) {
        if (cur_edge->info().geom_id == e_id && cur_edge->info().geom_dim == e_dim)
            classified_edges++;
    }
    for (auto cur_node : bl.get_all_nodes()) {
        if (cur_node->info().geom_id == e_id && cur_node->info().geom_dim == e_dim)
            classified_nodes++;
    }
    REQUIRE(classified_edges == 2);
    REQUIRE(classified_nodes == 1);

}

TEST_CASE("BlockingTestSuite - split_one_block_twice", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::cblock::Blocking bl(&geom_model, true);
    gecko::cblock::BlockingClassifier cl(&bl);
    cl.classify();

    auto e = bl.get_all_edges()[0];
    auto e_id = e->info().geom_id;
    auto e_dim = e->info().geom_dim;

    auto e2 = bl.cmap()->attribute<1>(bl.cmap()->beta<1>(e->dart()));
    bl.cut_sheet(e);
    REQUIRE(bl.get_nb_cells<0>() == 12);
    REQUIRE(bl.get_nb_cells<1>() == 20);
    REQUIRE(bl.get_nb_cells<2>() == 11);
    REQUIRE(bl.get_nb_cells<3>() == 2);
    REQUIRE(bl.cmap()->is_valid());

    int classified_nodes = 0;
    int classified_edges = 0;
    for (auto cur_edge : bl.get_all_edges()) {
        if (cur_edge->info().geom_id == e_id && cur_edge->info().geom_dim == e_dim)
            classified_edges++;
    }
    for (auto cur_node : bl.get_all_nodes()) {
        if (cur_node->info().geom_id == e_id && cur_node->info().geom_dim == e_dim)
            classified_nodes++;
    }
    REQUIRE(classified_edges == 2);
    REQUIRE(classified_nodes == 1);

    gmds::math::Point p_cut(-5, 5, 2);
    bl.cut_sheet(e2, p_cut);
    for (auto n : bl.get_all_nodes()) {
        auto nz = n->info().point.Z();
        REQUIRE((Approx(nz).margin(1e-4) == 5 || Approx(nz).margin(1e-4) == -5 || Approx(nz).margin(1e-4) == 2));
    }
    REQUIRE(bl.get_nb_cells<0>() == 18);
    REQUIRE(bl.get_nb_cells<1>() == 33);
    REQUIRE(bl.get_nb_cells<2>() == 20);
    REQUIRE(bl.get_nb_cells<3>() == 4);
}

TEST_CASE("BlockingTestSuite - split_until", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::cblock::Blocking bl(&geom_model, true);
    gecko::cblock::BlockingClassifier cl(&bl);
    cl.classify();

    bool to_cut = true;
    gecko::cblock::Blocking::Edge edge_to_cut;
    while (to_cut) {

        to_cut = false;
        auto all_edges = bl.get_all_edges();
        for (auto cur_edge : all_edges) {
            auto nodes_of_e = bl.get_nodes_of_edge(cur_edge);
            gmds::math::Point p0 = nodes_of_e[0]->info().point;
            gmds::math::Point p1 = nodes_of_e[1]->info().point;
            if (p0.distance(p1) > 2) {
                to_cut = true;
                edge_to_cut = cur_edge;
            }
        }
        if (to_cut) {
            bl.cut_sheet(edge_to_cut, 0.5);
        }
    }
    REQUIRE(bl.get_nb_cells<3>() == 512);
}

TEST_CASE("BlockingTestSuite - init_from_geom_bounding_box", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::cblock::Blocking bl(&geom_model, true);

    REQUIRE(bl.get_nb_cells<0>() == 8);
    REQUIRE(bl.get_nb_cells<1>() == 12);
    REQUIRE(bl.get_nb_cells<2>() == 6);
    REQUIRE(bl.get_nb_cells<3>() == 1);

    for (auto a : bl.cmap()->attributes<0>()) {
        gmds::math::Point p = a.info().point;
        REQUIRE(Approx(fabs(p.X())).margin(1e-8) == 5);
        REQUIRE(Approx(fabs(p.Y())).margin(1e-8) == 5);
        REQUIRE(Approx(fabs(p.Z())).margin(1e-8) == 5);
    }
}

TEST_CASE("BlockingTestSuite - single_block_to_mesh", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::cblock::Blocking bl(&geom_model);

    gmds::math::Point p000(0, 0, 0);
    gmds::math::Point p010(0, 1, 0);
    gmds::math::Point p110(1, 1, 0);
    gmds::math::Point p100(1, 0, 0);
    gmds::math::Point p001(0, 0, 1);
    gmds::math::Point p011(0, 1, 1);
    gmds::math::Point p111(1, 1, 1);
    gmds::math::Point p101(1, 0, 1);

    bl.create_block(p000, p010, p110, p100, p001, p011, p111, p101);

    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F | gmds::R | gmds::E2N | gmds::F2N | gmds::R2N));
    bl.convert_to_mesh(m);

    REQUIRE(m.getNbNodes() == 8);
    REQUIRE(m.getNbEdges() == 12);
    REQUIRE(m.getNbFaces() == 6);
    REQUIRE(m.getNbRegions() == 1);
}

TEST_CASE("BlockingTestSuite - test_topological_queries", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::cblock::Blocking bl(&geom_model, true);

    auto bl_nodes = bl.get_all_nodes();
    auto bl_edges = bl.get_all_edges();
    auto bl_faces = bl.get_all_faces();

    for (auto n : bl_nodes) {
        auto fs = bl.get_faces_of_node(n);
        REQUIRE(fs.size() == 3);

        auto es = bl.get_edges_of_node(n);
        REQUIRE(es.size() == 3);

        auto bs = bl.get_blocks_of_node(n);
        REQUIRE(bs.size() == 1);
    }

    for (auto e : bl_edges) {
        auto ns = bl.get_nodes_of_edge(e);
        REQUIRE(ns.size() == 2);

        auto fs = bl.get_faces_of_edge(e);
        REQUIRE(fs.size() == 2);

        auto bs = bl.get_blocks_of_edge(e);
        REQUIRE(bs.size() == 1);
    }

    for (auto f : bl_faces) {
        auto ns = bl.get_nodes_of_face(f);
        REQUIRE(ns.size() == 4);

        auto es = bl.get_edges_of_face(f);
        REQUIRE(es.size() == 4);

        auto bs = bl.get_blocks_of_face(f);
        REQUIRE(bs.size() == 1);
    }
}


TEST_CASE("BlockingTestSuite - test_init_from_ig_mesh", "[BlockingTestSuite]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::cblock::Blocking bl(&geom_model, false);

    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::N |  gmds::R | gmds::R2N));
    gmds::Node n0 = m.newNode(gmds::math::Point(0,0,0));
    gmds::Node n1 = m.newNode(gmds::math::Point(1,0,0));
    gmds::Node n2 = m.newNode(gmds::math::Point(1,1,0));
    gmds::Node n3 = m.newNode(gmds::math::Point(0,1,0));

    gmds::Node n4 = m.newNode(gmds::math::Point(0,0,1));
    gmds::Node n5 = m.newNode(gmds::math::Point(1,0,1));
    gmds::Node n6 = m.newNode(gmds::math::Point(1,1,1));
    gmds::Node n7 = m.newNode(gmds::math::Point(0,1,1));

    gmds::Node n8 = m.newNode(gmds::math::Point(0,0,2));
    gmds::Node n9 = m.newNode(gmds::math::Point(1,0,2));
    gmds::Node n10= m.newNode(gmds::math::Point(1,1,2));
    gmds::Node n11= m.newNode(gmds::math::Point(0,1,2));

    gmds::Node n12= m.newNode(gmds::math::Point(0,0,3));
    gmds::Node n13= m.newNode(gmds::math::Point(1,0,3));
    gmds::Node n14= m.newNode(gmds::math::Point(1,1,3));
    gmds::Node n15= m.newNode(gmds::math::Point(0,1,3));

    m.newHex(n0,n1,n2,n3,n4,n5,n6,n7);
    m.newHex(n4,n5,n6,n7,n8,n9,n10,n11);
    m.newHex(n8,n9,n10,n11,n12,n13,n14,n15);
    bl.init_from_mesh(m);

    REQUIRE(bl.get_all_nodes().size() == 16);
    REQUIRE(bl.get_all_edges().size() == 28);
    REQUIRE(bl.get_all_faces().size() == 16);
    REQUIRE(bl.get_all_blocks().size() == 3);
}
