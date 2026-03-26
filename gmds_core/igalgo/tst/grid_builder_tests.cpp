#include <catch2/catch_test_macros.hpp>
#include <unit_test_config.h>
#include <gmds/ig/Mesh.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>

TEST_CASE("GridBuildOpClass: 2D grid generation", "[GridBuilder][2D]") {
    // Création du maillage 2D
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::F | gmds::F2N));
    gmds::GridBuilder gb(&m, 2);

    REQUIRE(gb.isValid()); // Vérification de la validité du GridBuilder

    // Exécution de la construction de la grille
    gb.execute(3, 1.0, 4, 1.0);

    // Vérification des résultats
    REQUIRE(m.getNbNodes() == 12); // Nombre attendu de nœuds
    REQUIRE(m.getNbFaces() == 6);  // Nombre attendu de faces
}

TEST_CASE("GridBuildOpClass: 3D grid generation", "[GridBuilder][3D]") {
    // Création du maillage 3D
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::R | gmds::R2N));
    gmds::GridBuilder gb(&m, 3);

    REQUIRE(gb.isValid()); // Vérification de la validité du GridBuilder

    // Exécution de la construction de la grille
    gb.execute(3, 1.0, 4, 1.0, 3, 2.0);

    // Vérification des résultats
    REQUIRE(m.getNbNodes() == 36);    // Nombre attendu de nœuds
    REQUIRE(m.getNbRegions() == 12); // Nombre attendu de régions
}
