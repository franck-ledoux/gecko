#include <gecko/gblock/Blocking.h>
#include <gecko/gblock/BlockingClassifier.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <unit_test_config.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

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

std::tuple<int,int,int,int> get_node_statistics(gecko::gblock::Blocking& ABlocking){
    auto nb_on_vertex=0;
    auto nb_on_curve=0;
    auto nb_on_surface=0;
    auto nb_in_volume=0;
    std::vector<gecko::gblock::Blocking::Node> all_nodes = ABlocking.get_all_nodes();
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

std::tuple<int,int,int> get_edge_statistics(gecko::gblock::Blocking& ABlocking){
    auto nb_on_curve=0;
    auto nb_on_surface=0;
    auto nb_in_volume=0;
    std::vector<gecko::gblock::Blocking::Edge> all_edges = ABlocking.get_all_edges();
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

std::tuple<int,int> get_face_statistics(gecko::gblock::Blocking& ABlocking){
    auto nb_on_surface=0;
    auto nb_in_volume=0;
    std::vector<gecko::gblock::Blocking::Face> all_faces = ABlocking.get_all_faces();
    for(auto f:all_faces){
        if(f->info().geom_dim==2)
            nb_on_surface++;
        else if(f->info().geom_dim==3)
            nb_in_volume++;
    }
    return std::make_tuple(nb_on_surface,nb_in_volume);
}

void export_vtk(gecko::gblock::Blocking& ABlocking, int AModel, const std::string& AFileName){
    gmds::Mesh m_out(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F | gmds::R | gmds::E2N | gmds::F2N | gmds::R2N));
    ABlocking.convert_to_mesh(m_out);
    gmds::IGMeshIOService ioService(&m_out);
    gmds::VTKWriter writer(&ioService);
    writer.setCellOptions(AModel);
    writer.setDataOptions(AModel);
    writer.write(AFileName);
}

TEST_CASE("global_cell_accessors", "[BlockingTestSuite]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model);

    gmds::math::Point p000(0, 0, 0);
    gmds::math::Point p010(0, 1, 0);
    gmds::math::Point p110(1, 1, 0);
    gmds::math::Point p100(1, 0, 0);

    gmds::math::Point p001(0, 0, 1);
    gmds::math::Point p011(0, 1, 1);
    gmds::math::Point p111(1, 1, 1);
    gmds::math::Point p101(1, 0, 1);

    gmds::math::Point p002(0, 0, 2);
    gmds::math::Point p012(0, 1, 2);
    gmds::math::Point p112(1, 1, 2);
    gmds::math::Point p102(1, 0, 2);

    auto b1 = bl.create_block(p000, p010, p110, p100, p001, p011, p111, p101);
    auto b2 = bl.create_block(p001, p011, p111, p101, p002, p012, p112, p102);

    REQUIRE(bl.get_nb_cells<0>() == 16);
    REQUIRE(bl.get_nb_cells<1>() == 24);
    REQUIRE(bl.get_nb_cells<2>() == 12);
    REQUIRE(bl.get_nb_cells<3>() == 2);
    REQUIRE(bl.get_all_nodes().size() == 16);
    REQUIRE(bl.get_all_edges().size() == 24);
    REQUIRE(bl.get_all_faces().size() == 12);
    REQUIRE(bl.get_all_blocks().size() == 2);

    bl.sew<3>(b1->dart(), b2->dart());

    REQUIRE(bl.get_nb_cells<0>() == 12);
    REQUIRE(bl.get_nb_cells<1>() == 20);
    REQUIRE(bl.get_nb_cells<2>() == 11);
    REQUIRE(bl.get_nb_cells<3>() == 2);
    REQUIRE(bl.get_all_nodes().size() == 12);
    REQUIRE(bl.get_all_edges().size() == 20);
    REQUIRE(bl.get_all_faces().size() == 11);
    REQUIRE(bl.get_all_blocks().size() == 2);
}

TEST_CASE("remove_block", "[BlockingTestSuite]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model);

    gmds::math::Point p000(0, 0, 0);
    gmds::math::Point p010(0, 1, 0);
    gmds::math::Point p110(1, 1, 0);
    gmds::math::Point p100(1, 0, 0);

    gmds::math::Point p001(0, 0, 1);
    gmds::math::Point p011(0, 1, 1);
    gmds::math::Point p111(1, 1, 1);
    gmds::math::Point p101(1, 0, 1);

    gmds::math::Point p002(0, 0, 2);
    gmds::math::Point p012(0, 1, 2);
    gmds::math::Point p112(1, 1, 2);
    gmds::math::Point p102(1, 0, 2);

    auto b1 = bl.create_block(p000, p010, p110, p100, p001, p011, p111, p101);
    auto b2 = bl.create_block(p001, p011, p111, p101, p002, p012, p112, p102);
    bl.sew<3>(b1->dart(), b2->dart());

    REQUIRE(bl.get_nb_cells<3>() == 2);
    bl.remove_block(b1);
    REQUIRE(bl.get_nb_cells<3>() == 1);
}

TEST_CASE("single_block", "[BlockingTestSuite]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model);

    gmds::math::Point p000(0, 0, 0);
    gmds::math::Point p010(0, 1, 0);
    gmds::math::Point p110(1, 1, 0);
    gmds::math::Point p100(1, 0, 0);

    gmds::math::Point p001(0, 0, 1);
    gmds::math::Point p011(0, 1, 1);
    gmds::math::Point p111(1, 1, 1);
    gmds::math::Point p101(1, 0, 1);

    auto b = bl.create_block(p000, p010, p110, p100, p001, p011, p111, p101);

    REQUIRE(b->info().geom_dim == 3);
    REQUIRE(b->info().geom_id == gmds::NullID);

    auto fs = bl.get_faces_of_block(b);
    REQUIRE(fs.size() == 6);

    auto block_center = bl.get_center_of_block(b);
    for (auto f : fs) {
        gmds::math::Point face_center = bl.get_center_of_face(f);
        REQUIRE(std::abs(block_center.distance(face_center) - 0.5) < 1e-8);
    }
}

TEST_CASE("single_block_parallel_edges", "[BlockingTestSuite]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model);

    gmds::math::Point p000(0, 0, 0);
    gmds::math::Point p010(0, 1, 0);
    gmds::math::Point p110(1, 1, 0);
    gmds::math::Point p100(1, 0, 0);

    gmds::math::Point p001(0, 0, 1);
    gmds::math::Point p011(0, 1, 1);
    gmds::math::Point p111(1, 1, 1);
    gmds::math::Point p101(1, 0, 1);

    bl.create_block(p000, p010, p110, p100, p001, p011, p111, p101);

    for (auto it = bl.gmap()->attributes<1>().begin(), itend = bl.gmap()->attributes<1>().end(); it != itend; ++it) {
        std::vector<gecko::gblock::Blocking::Edge> parallel_edges;
        bl.get_all_sheet_edges(it, parallel_edges);
        REQUIRE(parallel_edges.size() == 4);
    }

    std::vector<std::vector<gecko::gblock::Blocking::Edge>> all_edges = bl.get_all_sheet_edge_sets();
    REQUIRE(all_edges.size() == 3);
    for (auto& sh_edges : all_edges) {
        REQUIRE(sh_edges.size() == 4);
    }
}

TEST_CASE("get_edges_of_a_block", "[BlockingTestSuite]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, true);

    auto b = bl.get_all_blocks()[0];
    auto edges = bl.get_edges_of_block(b);
    REQUIRE(edges.size() == 12);
}


using Catch::Approx;

TEST_CASE("BlockingTestSuite - split_one_block_twice", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, true);
    gecko::gblock::BlockingClassifier cl(&bl);
    cl.classify();

    auto e = bl.get_all_edges()[0];
    auto e_id = e->info().geom_id;
    auto e_dim = e->info().geom_dim;

    auto e2 = bl.gmap()->attribute<1>(bl.gmap()->alpha<1>(e->dart()));
    bl.cut_sheet(e);
    REQUIRE(bl.get_nb_cells<0>() == 12);
    REQUIRE(bl.get_nb_cells<1>() == 20);
    REQUIRE(bl.get_nb_cells<2>() == 11);
    REQUIRE(bl.get_nb_cells<3>() == 2);
    REQUIRE(bl.gmap()->is_valid());

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

TEST_CASE("BlockingTestSuite - cut_sheet_param_propag", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, true);
    gecko::gblock::BlockingClassifier cl(&bl);
    cl.classify();

    auto all_faces = bl.get_all_faces();
    std::vector<gecko::gblock::Blocking::Face> surf;
    for (auto f : all_faces) {
        gmds::math::Point ci = bl.get_center_of_face(f);
        if (fabs(ci.X() - 5) < 0.1 || fabs(ci.Y() - 5) < 0.1 || fabs(ci.Z() - 5) < 0.1) {
            surf.push_back(f);
        }
    }
    REQUIRE(bl.pillow(surf));

    bl.smooth(10);

    auto e = bl.get_edge(0, 3);
    bl.cut_sheet(e, 0.25);
    REQUIRE(bl.get_nb_cells<3>() == 7);
}

TEST_CASE("BlockingTestSuite - init_from_geom_bounding_box", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, true);

    REQUIRE(bl.get_nb_cells<0>() == 8);
    REQUIRE(bl.get_nb_cells<1>() == 12);
    REQUIRE(bl.get_nb_cells<2>() == 6);
    REQUIRE(bl.get_nb_cells<3>() == 1);

    for (auto a : bl.gmap()->attributes<0>()) {
        gmds::math::Point p = a.info().point;
        REQUIRE(Approx(fabs(p.X())).margin(1e-8) == 5);
        REQUIRE(Approx(fabs(p.Y())).margin(1e-8) == 5);
        REQUIRE(Approx(fabs(p.Z())).margin(1e-8) == 5);
    }
}

TEST_CASE("BlockingTestSuite - single_block_to_mesh", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model);

    gmds::math::Point p000(0, 0, 0);
    gmds::math::Point p010(0, 1, 0);
    gmds::math::Point p110(1, 1, 0);
    gmds::math::Point p100(1, 0, 0);
    gmds::math::Point p001(0, 0, 1);
    gmds::math::Point p011(0, 1, 1);
    gmds::math::Point p111(1, 1, 1);
    gmds::math::Point p101(1, 0, 1);

    auto b = bl.create_block(p000, p010, p110, p100, p001, p011, p111, p101);

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
    gecko::gblock::Blocking bl(&geom_model, true);

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
    gecko::gblock::Blocking bl(&geom_model, false);

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

TEST_CASE("BlockingTestSuite - test_chord_query", "[BlockingTestSuite]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, false);

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

    std::vector<gecko::gblock::Dart3> darts;
    std::vector<gecko::gblock::Blocking::Face> bl_faces = bl.get_all_faces();
    for (auto f : bl_faces) {
        gmds::math::Point center = bl.get_center_of_face(f);
        double z_integer_part = std::floor(center.Z());
        double z_decimal_part = center.Z() - z_integer_part;
        if (z_decimal_part == 0) {
            bl.get_all_chord_darts(f, darts);
            REQUIRE(darts.size() == 4);
            REQUIRE(bl.get_all_chord_blocks(f).size() == 3);
        } else {
            bl.get_all_chord_darts(f, darts);
            REQUIRE(darts.size() == 2);
            REQUIRE(bl.get_all_chord_blocks(f).size() == 1);
        }
    }
}

TEST_CASE("BlockingTestSuite - test_chord_collapse", "[BlockingTestSuite]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, false);

    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::R|gmds::R2N));
    gmds::GridBuilder gb(&m, 3);
    gb.execute(4, 1.0, 4, 1.0, 4, 1.0);
    bl.init_from_mesh(m);

    std::vector<gecko::gblock::Blocking::Face> bl_faces = bl.get_all_faces();
    REQUIRE(bl.get_nb_cells<3>() == 27);

    bool found_face = false;
    auto face_id = -1;
    gmds::math::Point seed(1.5, 1.5, 0.0);
    for (auto i = 0; i < bl_faces.size() && !found_face; i++) {
        gecko::gblock::Blocking::Face fi = bl_faces[i];
        gmds::math::Point ci = bl.get_center_of_face(fi);
        if (ci.distance2(seed) < 0.1) {
            found_face = true;
            face_id = i;
        }
    }

    std::vector<gecko::gblock::Blocking::Node> f_nodes = bl.get_nodes_of_face(bl_faces[face_id]);
    bl.collapse_chord(bl_faces[face_id], f_nodes[0], f_nodes[2]);

    REQUIRE(bl.get_nb_cells<3>() == 24);

    gmds::Mesh m_out(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F | gmds::R | gmds::E2N | gmds::F2N | gmds::R2N));
    bl.convert_to_mesh(m_out);

    gmds::IGMeshIOService ioService(&m_out);
    gmds::VTKWriter writer(&ioService);
    writer.setCellOptions(gmds::N | gmds::F);
    writer.setDataOptions(gmds::N | gmds::F);
    writer.write("collapse.vtk");
}


TEST_CASE("BlockingTestSuite test_pillow_6", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, true);
    gecko::gblock::BlockingClassifier cl(&bl);
    cl.classify();

    std::vector<gecko::gblock::Blocking::Face> all_faces = bl.get_all_faces();
    std::vector<gecko::gblock::Blocking::Face> surf;

    surf.clear();
    for (auto f : all_faces) {
        gmds::math::Point ci = bl.get_center_of_face(f);
        if (fabs(ci.X() - 5) < 0.1) {
            surf.push_back(f);
        } else if (fabs(ci.Y() - 5) < 0.1) {
            surf.push_back(f);
        } else if (fabs(ci.Y() + 5) < 0.1) {
            surf.push_back(f);
        } else if (fabs(ci.Z() - 5) < 0.1) {
            surf.push_back(f);
        } else if (fabs(ci.Z() + 5) < 0.1) {
            surf.push_back(f);
        }
    }

    REQUIRE(bl.pillow(surf));

    bl.smooth(10);

    REQUIRE(bl.get_nb_cells<3>() == 6);
    auto [nb_nodes_on_vertex, nb_nodes_on_curve, nb_nodes_on_surface, nb_nodes_in_volume] = get_node_statistics(bl);
    REQUIRE(nb_nodes_on_vertex == 8);
    REQUIRE(nb_nodes_on_curve == 0);
    REQUIRE(nb_nodes_on_surface == 4);
    REQUIRE(nb_nodes_in_volume == 4);
    auto [nb_edges_on_curve, nb_edges_on_surface, nb_edges_in_volume] = get_edge_statistics(bl);
    REQUIRE(nb_edges_on_curve == 12);
    REQUIRE(nb_edges_on_surface == 8);
    REQUIRE(nb_edges_in_volume == 12);
    auto [nb_faces_on_surface, nb_faces_in_volume] = get_face_statistics(bl);
    REQUIRE(nb_faces_on_surface == 10);
    REQUIRE(nb_faces_in_volume == 13);
}

TEST_CASE("BlockingTestSuite test_pillow_7", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, true);
    gecko::gblock::BlockingClassifier cl(&bl);
    cl.classify();

    std::vector<gecko::gblock::Blocking::Face> all_faces = bl.get_all_faces();

    REQUIRE(bl.pillow(all_faces));

    bl.smooth(10);

    REQUIRE(bl.get_nb_cells<3>() == 7);
    auto [nb_nodes_on_vertex, nb_nodes_on_curve, nb_nodes_on_surface, nb_nodes_in_volume] = get_node_statistics(bl);
    REQUIRE(nb_nodes_on_vertex == 8);
    REQUIRE(nb_nodes_on_curve == 0);
    REQUIRE(nb_nodes_on_surface == 0);
    REQUIRE(nb_nodes_in_volume == 8);
    auto [nb_edges_on_curve, nb_edges_on_surface, nb_edges_in_volume] = get_edge_statistics(bl);
    REQUIRE(nb_edges_on_curve == 12);
    REQUIRE(nb_edges_on_surface == 0);
    REQUIRE(nb_edges_in_volume == 20);
    auto [nb_faces_on_surface, nb_faces_in_volume] = get_face_statistics(bl);
    REQUIRE(nb_faces_on_surface == 6);
    REQUIRE(nb_faces_in_volume == 18);
}

TEST_CASE("BlockingTestSuite test_pillow_8", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, true);
    gecko::gblock::BlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gecko::gblock::Blocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for (auto par_edges_i : edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    REQUIRE(bl.get_nb_cells<3>() == 8);

    std::vector<gecko::gblock::Blocking::Face> all_faces = bl.get_all_faces();
    std::vector<gecko::gblock::Blocking::Face> surf;

    surf.clear();
    for (auto f : all_faces) {
        auto nb_adj = bl.get_blocks_of_face(f);
        if (nb_adj.size() == 1) {
            surf.push_back(f);
        }
    }

    REQUIRE(bl.pillow(surf));

    bl.smooth(10);

    REQUIRE(bl.get_nb_cells<3>() == 32);
    auto [nb_nodes_on_vertex, nb_nodes_on_curve, nb_nodes_on_surface, nb_nodes_in_volume] = get_node_statistics(bl);
    REQUIRE(nb_nodes_on_vertex == 8);
    REQUIRE(nb_nodes_on_curve == 12);
    REQUIRE(nb_nodes_on_surface == 6);
    REQUIRE(nb_nodes_in_volume == 27);
    auto [nb_edges_on_curve, nb_edges_on_surface, nb_edges_in_volume] = get_edge_statistics(bl);
    REQUIRE(nb_edges_on_curve == 24);
    REQUIRE(nb_edges_on_surface == 24);
    REQUIRE(nb_edges_in_volume == 80);
    auto [nb_faces_on_surface, nb_faces_in_volume] = get_face_statistics(bl);
    REQUIRE(nb_faces_on_surface == 24);
    REQUIRE(nb_faces_in_volume == 84);
}

TEST_CASE("BlockingTestSuite test_pillow_9", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, true);
    gecko::gblock::BlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gecko::gblock::Blocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for (auto par_edges_i : edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    REQUIRE(bl.get_nb_cells<3>() == 8);

    std::vector<gecko::gblock::Blocking::Face> all_faces = bl.get_all_faces();
    std::vector<gecko::gblock::Blocking::Face> surf;

    surf.clear();
    for (auto f : all_faces) {
        gmds::math::Point ci = bl.get_center_of_face(f);
        if (fabs(ci.X() - 5) < 0.1) {
            surf.push_back(f);
        } else if (fabs(ci.Y() - 5) < 0.1) {
            surf.push_back(f);
        } else if (fabs(ci.Z() - 5) < 0.1) {
            surf.push_back(f);
        }
    }

    REQUIRE(bl.pillow(surf));

    bl.smooth(10);

    REQUIRE(bl.get_nb_cells<3>() == 20);
    auto [nb_nodes_on_vertex, nb_nodes_on_curve, nb_nodes_on_surface, nb_nodes_in_volume] = get_node_statistics(bl);
    REQUIRE(nb_nodes_on_vertex == 8);
    REQUIRE(nb_nodes_on_curve == 15);
    REQUIRE(nb_nodes_on_surface == 15);
    REQUIRE(nb_nodes_in_volume == 8);
    auto [nb_edges_on_curve, nb_edges_on_surface, nb_edges_in_volume] = get_edge_statistics(bl);
    REQUIRE(nb_edges_on_curve == 27);
    REQUIRE(nb_edges_on_surface == 45);
    REQUIRE(nb_edges_in_volume == 31);
    auto [nb_faces_on_surface, nb_faces_in_volume] = get_face_statistics(bl);
    REQUIRE(nb_faces_on_surface == 36);
    REQUIRE(nb_faces_in_volume == 42);
}

TEST_CASE("BlockingTestSuite test_pillow_10", "[blocking]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, true);
    gecko::gblock::BlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gecko::gblock::Blocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for (auto par_edges_i : edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    REQUIRE(bl.get_nb_cells<3>() == 8);

    std::vector<gecko::gblock::Blocking::Face> all_faces = bl.get_all_faces();
    std::vector<gecko::gblock::Blocking::Face> surf;

    surf.clear();
    for (auto f : all_faces) {
        gmds::math::Point ci = bl.get_center_of_face(f);
        if (fabs(ci.X()) < 0.1) {
            surf.push_back(f);
        }
    }
    REQUIRE(surf.size() == 4);

    REQUIRE(bl.pillow(surf));

    bl.smooth(10);

    export_vtk(bl, gmds::N | gmds::E, "pillow_10_surf.vtk");

    REQUIRE(bl.get_nb_cells<3>() == 12);
    auto [nb_nodes_on_vertex, nb_nodes_on_curve, nb_nodes_on_surface, nb_nodes_in_volume] = get_node_statistics(bl);
    REQUIRE(nb_nodes_on_vertex == 8);
    REQUIRE(nb_nodes_on_curve == 16);
    REQUIRE(nb_nodes_on_surface == 10);
    REQUIRE(nb_nodes_in_volume == 2);
    auto [nb_edges_on_curve, nb_edges_on_surface, nb_edges_in_volume] = get_edge_statistics(bl);
    REQUIRE(nb_edges_on_curve == 28);
    REQUIRE(nb_edges_on_surface == 36);
    REQUIRE(nb_edges_in_volume == 11);
    auto [nb_faces_on_surface, nb_faces_in_volume] = get_face_statistics(bl);
    REQUIRE(nb_faces_on_surface == 32);
    REQUIRE(nb_faces_in_volume == 20);
}

TEST_CASE("test_pillow_11", "[BlockingTestSuite]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, true);
    gecko::gblock::BlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gecko::gblock::Blocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for (auto par_edges_i : edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    REQUIRE(bl.get_nb_cells<3>() == 8);

    std::vector<gecko::gblock::Blocking::Face> all_faces = bl.get_all_faces();
    std::vector<gecko::gblock::Blocking::Face> surf;
    surf.clear();

    for (auto f : all_faces) {
        gmds::math::Point ci = bl.get_center_of_face(f);
        if (fabs(ci.X()) < 0.1 && ci.Y() < 0 && ci.Z() < 0) {
            surf.push_back(f);
        } else if (ci.X() < 0 && fabs(ci.Y()) < 0.1 && ci.Z() < 0) {
            surf.push_back(f);
        } else if (ci.X() < 0 && ci.Y() < 0 && fabs(ci.Z()) < 0.1) {
            surf.push_back(f);
        }
    }
    REQUIRE(surf.size() == 3);

    REQUIRE(bl.pillow(surf));

    bl.smooth(10);

    REQUIRE(bl.get_nb_cells<3>() == 11);

    auto [nb_nodes_on_vertex, nb_nodes_on_curve, nb_nodes_on_surface, nb_nodes_in_volume] = get_node_statistics(bl);
    REQUIRE(nb_nodes_on_vertex == 8);
    REQUIRE(nb_nodes_on_curve == 15);
    REQUIRE(nb_nodes_on_surface == 9);
    REQUIRE(nb_nodes_in_volume == 2);

    auto [nb_edges_on_curve, nb_edges_on_surface, nb_edges_in_volume] = get_edge_statistics(bl);
    REQUIRE(nb_edges_on_curve == 27);
    REQUIRE(nb_edges_on_surface == 33);
    REQUIRE(nb_edges_in_volume == 10);

    auto [nb_faces_on_surface, nb_faces_in_volume] = get_face_statistics(bl);
    REQUIRE(nb_faces_on_surface == 30);
    REQUIRE(nb_faces_in_volume == 18);
}

TEST_CASE("test_pillow_12", "[BlockingTestSuite]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, true);
    gecko::gblock::BlockingClassifier cl(&bl);
    cl.classify();

    std::vector<std::vector<gecko::gblock::Blocking::Edge>> edges = bl.get_all_sheet_edge_sets();
    for (auto par_edges_i : edges) {
        bl.cut_sheet(par_edges_i[0]);
    }

    REQUIRE(bl.get_nb_cells<3>() == 8);

    std::vector<gecko::gblock::Blocking::Face> all_faces = bl.get_all_faces();
    std::vector<gecko::gblock::Blocking::Face> surf;
    surf.clear();

    for (auto f : all_faces) {
        gmds::math::Point ci = bl.get_center_of_face(f);
        if (fabs(ci.X()) < 0.1) {
            surf.push_back(f);
        } else if (ci.X() < 0 && fabs(ci.Y() - 5) < 0.1) {
            surf.push_back(f);
        }
    }
    REQUIRE(surf.size() == 6);

    REQUIRE(bl.pillow(surf));

    bl.smooth(10);

    REQUIRE(bl.get_nb_cells<3>() == 14);

    auto [nb_nodes_on_vertex, nb_nodes_on_curve, nb_nodes_on_surface, nb_nodes_in_volume] = get_node_statistics(bl);
    REQUIRE(nb_nodes_on_vertex == 8);
    REQUIRE(nb_nodes_on_curve == 16);
    REQUIRE(nb_nodes_on_surface == 12);
    REQUIRE(nb_nodes_in_volume == 3);

    auto [nb_edges_on_curve, nb_edges_on_surface, nb_edges_in_volume] = get_edge_statistics(bl);
    REQUIRE(nb_edges_on_curve == 28);
    REQUIRE(nb_edges_on_surface == 40);
    REQUIRE(nb_edges_in_volume == 15);

    auto [nb_faces_on_surface, nb_faces_in_volume] = get_face_statistics(bl);
    REQUIRE(nb_faces_on_surface == 34);
    REQUIRE(nb_faces_in_volume == 25);
}

TEST_CASE("save_vtk_blocking", "[BlockingTestSuite]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, true);
    gecko::gblock::BlockingClassifier classifier(&bl);

    bl.save_vtk_blocking("testSaveWork.vtk");
}

TEST_CASE("get_Id_block", "[BlockingTestSuite]") {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, true);
    gecko::gblock::BlockingClassifier classifier(&bl);

    classifier.clear_classification();
    auto errors = classifier.classify();

    auto listBlocks = bl.get_all_blocks();
    auto block = bl.get_block(listBlocks[0]->info().topo_id);

}
