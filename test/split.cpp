
#include <gecko/gblock/Blocking.h>
#include <gecko/gblock/BlockingClassifier.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>

#include <unit_test_config.h>

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

void export_vtk(gecko::gblock::Blocking& ABlocking, int AModel, const std::string& AFileName)
{
    gmds::Mesh m_out(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F | gmds::R | gmds::E2N | gmds::F2N | gmds::R2N));
    ABlocking.convert_to_mesh(m_out);
    gmds::IGMeshIOService ioService(&m_out);
    gmds::VTKWriter writer(&ioService);
    writer.setCellOptions(AModel);
    writer.setDataOptions(AModel);
    writer.write(AFileName);
}

int main() {
    gmds::cad::FACManager geom_model;
    setUp(geom_model);
    gecko::gblock::Blocking bl(&geom_model, true);
    gecko::gblock::BlockingClassifier cl(&bl);
    cl.classify();

    bool to_cut = true;
    gecko::gblock::Blocking::Edge edge_to_cut;
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

    export_vtk(bl,N|R, "multiple_cuts.vtk");
}