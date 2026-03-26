/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/igalgo/BoundaryExtractor3D.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cadfac/FACSurface.h>
#include <gmds/cadfac/FACVolume.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/IGMeshIOService.h>
#include <unit_test_config.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
/*----------------------------------------------------------------------------*/
using namespace gmds;

TEST_CASE("FACManager fromSurfMesh", "[FACManager]") {
    gmds::Mesh m_vol(gmds::MeshModel(DIM3 | R | F | E | N |
                                     R2N | R2F | R2E |
                                     F2N | F2R | F2E |
                                     E2F | E2N | N2E));
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

    cad::FACManager manager;
    manager.initFrom3DMesh(&m_vol);

    REQUIRE(manager.getNbPoints() == 8);
    REQUIRE(manager.getNbSurfaces() == 6);
    REQUIRE(manager.getNbCurves() == 12);

    std::vector<cad::GeomSurface*> surfs;
    manager.getSurfaces(surfs);

    for (auto s : surfs) {
        REQUIRE(s->computeArea() == Catch::Approx(100).margin(0.1));
        math::Point c(0, 0, 0);
        math::Point p = s->closestPoint(c);
        REQUIRE(c.distance(s->closestPoint(c)) == 5);
        REQUIRE(s->loops().size() == 1);
    }

    std::vector<cad::GeomCurve*> curvs;
    manager.getCurves(curvs);

    for (auto c : curvs) {
        REQUIRE(c->length() == 10);
        math::Point p(0, 0, 0);
        REQUIRE(p.distance(c->closestPoint(p)) == Catch::Approx(7.07107).margin(0.01));

        math::Point p_far(100000000, 0, 0);
        math::Point p_on = c->closestPoint(p_far);
        REQUIRE(p_on.distance(math::Point(0, 0, 0)) < 10);

        c->project(p_far);
        REQUIRE(p_far.distance(math::Point(0, 0, 0)) < 10);
        REQUIRE(p_on.distance(p_far) == Catch::Approx(0).margin(1e-9));
    }

    std::vector<cad::GeomPoint*> pnts;
    manager.getPoints(pnts);
    for (auto p : pnts) {
        math::Point c(0, 0, 0);
        REQUIRE(c.distance(p->point()) == Catch::Approx(8.66).margin(0.01));
    }

    REQUIRE(manager.getPoints().size() == 8);
    REQUIRE(manager.getCurves().size() == 12);
    REQUIRE(manager.getSurfaces().size() == 6);
    REQUIRE(manager.getVolumes().size() == 1);
    cad::GeomVolume* vol = manager.getVolumes()[0];
    REQUIRE(vol->isIn(math::Point(0,0,0)));
}

TEST_CASE("FACManager surf_projection", "[FACManager]") {
    gmds::Mesh m_vol(gmds::MeshModel(DIM3 | R | F | E | N |
                                     R2N | R2F | R2E |
                                     F2N | F2R | F2E |
                                     E2F | E2N | N2E));
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

    cad::FACManager manager;
    manager.initFrom3DMesh(&m_vol);

    cad::GeomSurface* s1 = manager.getSurface(1);

    TCoord min_s1[3], max_s1[3];
    s1->computeBoundingBox(min_s1, max_s1);

    math::Point p(-10.6674, 24.1967, 109.342);
    s1->project(p);
}

TEST_CASE("FACManager project", "[FACManager]") {
    gmds::Mesh m_vol(gmds::MeshModel(DIM3 | R | F | E | N |
                                     R2N | R2F | R2E |
                                     F2N | F2R | F2E |
                                     E2F | E2N | N2E));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir + "/simple_cube.vtk";

    gmds::IGMeshIOService ioService(&m_vol);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N | gmds::R);
    vtkReader.read(vtk_file);

    gmds::MeshDoctor doc(&m_vol);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    cad::FACManager manager;
    manager.initFrom3DMesh(&m_vol);

    cad::GeomSurface* s = manager.getSurface(2);

    TCoord min_s[3], max_s[3];
    s->computeBoundingBox(min_s, max_s);

    math::Point ps(0.1, -0.1, 0.1);
    math::Point ps_on(0.1, 0., 0.1);

    s->project(ps);
    REQUIRE(ps.distance(ps_on) == Catch::Approx(0).margin(1e-15));

    cad::GeomCurve* c = manager.getCurve(1);

    TCoord min_c[3], max_c[3];
    c->computeBoundingBox(min_c, max_c);

    math::Point pc(-0.1, 0.2, 3);
    math::Point pc_on(0., 0.2, 0.);

    c->project(pc);
    REQUIRE(pc.distance(pc_on) == Catch::Approx(0).margin(1e-15));
}

TEST_CASE("FACManager boundingbox", "[FACManager]") {
    gmds::Mesh m_vol(gmds::MeshModel(DIM3 | R | F | E | N |
                                     R2N | R2F | R2E |
                                     F2N | F2R | F2E |
                                     E2F | E2N | N2E));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir + "/B1_tet.vtk";

    gmds::IGMeshIOService ioService(&m_vol);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N | gmds::R);
    vtkReader.read(vtk_file);

    gmds::MeshDoctor doc(&m_vol);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    cad::FACManager manager;
    manager.initFrom3DMesh(&m_vol);

    cad::GeomVolume* v0 = manager.getVolume(1);
    REQUIRE(v0->id() == 1);

    TCoord minXYZ[3], maxXYZ[3];
    v0->computeBoundingBox(minXYZ, maxXYZ);
    REQUIRE(minXYZ[0] == Catch::Approx(-245.405));
    REQUIRE(minXYZ[1] == Catch::Approx(-245.407));
    REQUIRE(minXYZ[2] == Catch::Approx(-122.702));
    REQUIRE(maxXYZ[0] == Catch::Approx(638.053));
    REQUIRE(maxXYZ[1] == Catch::Approx(245.405));
    REQUIRE(maxXYZ[2] == Catch::Approx(368.107));
}

TEST_CASE("is in test", "[FACManager]") {
    gmds::Mesh vol(gmds::MeshModel(DIM3 | R | F | E | N |
                                     R2N | R2F | R2E |
                                     F2N | F2R | F2E |
                                     E2F | E2N | N2E));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir + "/cube_with_holes.vtk";

    gmds::IGMeshIOService ioService(&vol);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N | gmds::R);
    vtkReader.read(vtk_file);

    gmds::MeshDoctor doc(&vol);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    cad::FACManager manager;
    manager.initFrom3DMesh(&vol);

    cad::GeomVolume* gv = manager.getVolumes()[0];
    REQUIRE(gv->isIn(math::Point(4,4,4)) == false);
    REQUIRE(gv->isIn(math::Point(0,0,0)) == false);
    REQUIRE(gv->isIn(math::Point(4,0,0)) == true);
}