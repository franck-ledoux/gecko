/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
#include <unit_test_config.h>
#include <catch2/catch_test_macros.hpp>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST_CASE("from surface mesh", "[Geom Topology test]")
    {
        // WE WRITE
        Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                                   R2N|R2F|R2E|
                                   F2N|F2R|F2E|
                                   E2F|E2N|N2E|N2R));

        std::string dir(TEST_SAMPLES_DIR);
        std::string vtk_file = dir+"/tet_in_box.vtk";
        IGMeshIOService ioService(&m_vol);

        VTKReader vtkReader(&ioService);
        vtkReader.setCellOptions(gmds::N|gmds::R);
        vtkReader.read(vtk_file);

        MeshDoctor doc(&m_vol);
        doc.buildFacesAndR2F();
        doc.buildEdgesAndX2E();
        doc.updateUpwardConnectivity();


        cad::FACManager manager;
        cad::GeomMeshLinker linker;

        manager.initAndLinkFrom3DMesh(&m_vol,&linker);
        VTKWriter vtkW(&ioService);
        vtkW.setCellOptions(gmds::N|gmds::R);
        vtkW.setDataOptions(gmds::N|gmds::R);
        vtkW.write("toto_link.vtk");

        REQUIRE(cad::GeomMeshLinker::LinkPoint == linker.getGeomDim<Node>(1));
        REQUIRE(cad::GeomMeshLinker::LinkPoint == linker.getGeomDim<Node>(2));
        REQUIRE(cad::GeomMeshLinker::LinkPoint == linker.getGeomDim<Node>(4));
        REQUIRE(cad::GeomMeshLinker::LinkPoint == linker.getGeomDim<Node>(7));


        REQUIRE(cad::GeomMeshLinker::LinkCurve==linker.getGeomDim<Node>(9));
        REQUIRE(cad::GeomMeshLinker::LinkCurve== linker.getGeomDim<Node>(15));
        REQUIRE(cad::GeomMeshLinker::LinkCurve== linker.getGeomDim<Node>(18));
        REQUIRE(cad::GeomMeshLinker::LinkCurve== linker.getGeomDim<Node>(19));

        REQUIRE(cad::GeomMeshLinker::LinkSurface== linker.getGeomDim<Node>(20));
        REQUIRE(cad::GeomMeshLinker::LinkSurface== linker.getGeomDim<Node>(24));


        REQUIRE(2== linker.getGeomId<Node>(1));
        REQUIRE(3== linker.getGeomId<Node>(2));
        REQUIRE(5== linker.getGeomId<Node>(4));
        REQUIRE(8== linker.getGeomId<Node>(7));

        REQUIRE(4== linker.getGeomId<Node>(9));
        REQUIRE(9== linker.getGeomId<Node>(15));
        REQUIRE(11== linker.getGeomId<Node>(18));
        REQUIRE(2== linker.getGeomId<Node>(19));

        REQUIRE(4== linker.getGeomId<Node>(20));
        REQUIRE(3== linker.getGeomId<Node>(24));

        REQUIRE(cad::GeomMeshLinker::LinkPoint  == linker.getGeomInfo<Node>(4). first);
        REQUIRE(cad::GeomMeshLinker::LinkCurve  == linker.getGeomInfo<Node>(18).first);
        REQUIRE(cad::GeomMeshLinker::LinkSurface== linker.getGeomInfo<Node>(24).first);
        REQUIRE(5 == linker.getGeomInfo<Node>(4). second);
        REQUIRE(11== linker.getGeomInfo<Node>(18).second);
        REQUIRE(3 == linker.getGeomInfo<Node>(24).second);

        Node n1 = m_vol.get<Node>(1);
        REQUIRE(cad::GeomMeshLinker::LinkPoint==linker.getGeomDim(n1));
        REQUIRE(2== linker.getGeomId(n1));
    }
