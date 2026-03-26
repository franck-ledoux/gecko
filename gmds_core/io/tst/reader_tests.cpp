#include <catch2/catch_test_macros.hpp>

#include <gmds/ig/Mesh.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/MeditReader.h>
#include <gmds/io/IGMeshIOService.h>

#include <unit_test_config.h>


TEST_CASE("MeshClass: testIOVTK_NF", "[Mesh][VTK]") {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir + "/triangulated_quad.vtk";

    gmds::IGMeshIOService ioService(&m);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N | gmds::F);
    vtkReader.read(vtk_file);

    REQUIRE(m.getNbNodes() == 29);
    REQUIRE(m.getNbFaces() == 40);
}

TEST_CASE("MeshClass: testIOVTK_NR", "[Mesh][VTK]") {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::N));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir + "/tet_in_box.vtk";

    gmds::IGMeshIOService ioService(&m);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N | gmds::R);
    vtkReader.read(vtk_file);

    REQUIRE(m.getNbNodes() == 29);
    REQUIRE(m.getNbRegions() == 61);
}

TEST_CASE("MeshClass: testIOVTK_NEF", "[Mesh][VTK]") {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::E | gmds::N));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir + "/triangulated_quad.vtk";

    gmds::IGMeshIOService ioService(&m);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N | gmds::E | gmds::F);
    vtkReader.read(vtk_file);

    REQUIRE(m.getNbNodes() == 29);
    REQUIRE(m.getNbEdges() == 16);
    REQUIRE(m.getNbFaces() == 40);
}

TEST_CASE("MeshClass: testIOMedit", "[Mesh][Medit]") {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));

    std::string dir(TEST_SAMPLES_DIR);
    std::string med_file = dir + "/triangulated_quad.mesh";

    gmds::IGMeshIOService ioService(&m);
    gmds::MeditReader meditReader(&ioService);

    meditReader.setCellOptions(gmds::N | gmds::F);
    meditReader.read(med_file);

    REQUIRE(m.getNbNodes() == 29);
    REQUIRE(m.getNbFaces() == 40);
}

TEST_CASE("MeshClass: testIOVTK_several_types", "[Mesh][VTK]") {
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::N | gmds::R2N));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir + "/several_types.vtk";

    gmds::IGMeshIOService ioService(&m);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N | gmds::R);
    vtkReader.read(vtk_file);

    REQUIRE(m.getNbNodes() == 14);
    REQUIRE(m.getNbRegions() == 4);
    REQUIRE(m.getNbTetrahedra() == 1);
    REQUIRE(m.getNbHexahedra() == 2);
}
