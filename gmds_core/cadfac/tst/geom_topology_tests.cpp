/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/igalgo/GridBuilder.h>
#include <unit_test_config.h>
#include <gmds/io/IGMeshIOService.h>
#include <catch2/catch_test_macros.hpp>
#include <gmds/io/VTKReader.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST_CASE("connectivity", "[Geom Topology test]") {

  Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                               R2N|R2F|R2E|
                               F2N|F2R|F2E|
                               E2F|E2N|N2E|N2R|N2F));

    GridBuilder gb(&m_vol, 3);
    gb.execute(3,1.,3,1.,3,1.);
    REQUIRE(m_vol.getNbRegions()==8);
    MeshDoctor doc(&m_vol);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    cad::FACManager manager;
    cad::GeomMeshLinker linker;

    manager.initAndLinkFrom3DMesh(&m_vol,&linker);

//==================================
// tests on the points
    std::vector<cad::GeomPoint*> pnts;
    manager.getPoints(pnts);
    for(auto p:pnts){
        REQUIRE(p->volumes().size()==1);
        REQUIRE(p->surfaces().size()==3);
        REQUIRE(p->curves().size()==3);
    }
//==================================
// tests on the curves
    std::vector<cad::GeomCurve*> curves;
    manager.getCurves(curves);
    for(auto c:curves){
        REQUIRE(c->volumes().size()==1);
        std::vector<cad::GeomSurface*> surfs =c->surfaces();
        std::vector<cad::GeomPoint*> pts =c->points();
        REQUIRE(c->surfaces().size()==2);
        REQUIRE(c->points().size()==2);
    }
//==================================
// tests on the surfaces
    std::vector<cad::GeomSurface*> surfs;
    manager.getSurfaces(surfs);
    for(auto s:surfs){
        REQUIRE(s->volumes().size()==1);
        REQUIRE(s->curves().size()==4);
        REQUIRE(s->points().size()==4);
    }
//==================================
// tests on the volumes
    std::vector<cad::GeomVolume*> vols;
    manager.getVolumes(vols);
    for(auto v:vols){
        REQUIRE(v->surfaces().size()==6);
        REQUIRE(v->curves().size()==12);
        REQUIRE(v->points().size()==8);
    }
}
/*----------------------------------------------------------------------------*/
TEST_CASE("test_convex_curve", "[Geom Topology test]") {
    Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                               R2N|R2F|R2E|
                               F2N|F2R|F2E|
                               E2F|E2N|N2E|N2R|N2F));

    GridBuilder gb(&m_vol, 3);
    gb.execute(3,1.,3,1.,3,1.);
    REQUIRE(m_vol.getNbRegions()==8);
    MeshDoctor doc(&m_vol);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    cad::FACManager manager;
    cad::GeomMeshLinker linker;

    manager.initAndLinkFrom3DMesh(&m_vol,&linker);

    manager.write_surfaces("fac_model.vtk");
    // Check that all the cuvrve are strictly convex
    std::vector<cad::GeomCurve*> curves;
    manager.getCurves(curves);
    for(auto c:curves){
        REQUIRE(c->getCurvatureInfo()==  cad::GeomCurve::Convex);
    }
}
/*----------------------------------------------------------------------------*/
TEST_CASE("cube convex B45", "[Geom Topology test]") {
    Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                               R2N|R2F|R2E|
                               F2N|F2R|F2E|
                               E2F|E2N|N2E|N2R|N2F));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/B45.vtk";
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

    manager.write_surfaces("fac_model.vtk");
    // Check that all the cuvrve are strictly convex
    std::vector<cad::GeomCurve*> curves;
    manager.getCurves(curves);
    int nb_convex =0;
    int nb_concave=0;
    for(auto c:curves){
        cad::GeomCurve::CurvatureInfo ci = c->getCurvatureInfo();
        if(ci==cad::GeomCurve::Convex)
            nb_convex++;
        else if(ci==cad::GeomCurve::Concave)
            nb_concave++;
    }
    REQUIRE(17== nb_convex);
    REQUIRE(1 == nb_concave);
}
