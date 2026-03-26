#include <catch2/catch_test_macros.hpp>
#include <unit_test_config.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/igalgo/SurfaceReorient.h>
#include <iostream>
#include <gmds/io/IGMeshIOService.h>

using namespace gmds;

TEST_CASE("SurfaceReorientTestClass: testReorient", "[Surface][Reorient]") {
	// Initialisation de la grille
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N | gmds::N2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/reorient3D_test.vtk";

	// Chargement des données
	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	vtkReader.read(vtk_file);

	// Préparation du maillage
	gmds::MeshDoctor doc(&m);
	doc.updateUpwardConnectivity();

	// Reorientation des surfaces
	gmds::SurfaceReorient reorient(&m);
	REQUIRE(reorient.isValid());

	auto nb_reorient = reorient.execute();
	REQUIRE(nb_reorient != 0);

	// Écriture des résultats
	gmds::VTKWriter writer(&ioService);
	writer.setCellOptions(gmds::N | gmds::F);
	writer.write("toto.vtk");
}
