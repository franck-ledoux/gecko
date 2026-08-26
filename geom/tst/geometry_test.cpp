#include <algorithm>
#include <filesystem>
#include <ranges>
#include <type_traits>
#include <variant>

#include <unit_test_config.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <gecko/geom/FacetedGeometry.h>
#include <gecko/io/GmshMeshWriter.h>
#include <catch2/matchers/catch_matchers_range_equals.hpp>

using Catch::Approx;
using namespace gecko;

namespace {
    /** @brief Path used to write a fixture .msh file read back by a test. */
    std::filesystem::path fixture_path(const std::string &name) {
        return std::filesystem::temp_directory_path() / name;
    }

    /**
     * @brief Finds the entity in @p range whose entity_tag() matches @p tag.
     * @throw Catch2 assertion failure (via REQUIRE) if none is found.
     */
    template<typename Range>
    const auto &find_by_tag(const Range &range, Int tag) {
        auto it = std::ranges::find_if(range, [tag](const auto &e) { return e.entity_tag() == tag; });
        REQUIRE(it != std::ranges::end(range));
        return *it;
    }
} // namespace

TEST_CASE("FacetedGeometry_BoundaryRep_FromTriangles", "[FacetedGeometry]") {
    SimplicialMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(1, 1, 0);
    auto n3 = mesh.add_node(0, 1, 0);

    GroupRegistry groups;
    auto wall = groups.add_group("Wall", GroupDim::Dim2);
    auto edges_group = groups.add_group("Edges", GroupDim::Dim1);
    auto surfaces_group = groups.add_group("Surfaces", GroupDim::Dim2);

    auto &node_group = mesh.add_variable<GroupId, CellType::Node>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    auto &node_entity = mesh.add_variable<Int, CellType::Node>(std::string(io::ENTITY_TAG_VARIABLE));
    auto &edge_group_var = mesh.add_variable<GroupId, CellType::Edge>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    auto &edge_entity_var = mesh.add_variable<Int, CellType::Edge>(std::string(io::ENTITY_TAG_VARIABLE));
    auto &face_group_var = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    auto &face_entity_var = mesh.add_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));

    // Only n0 is an explicit (tagged) Gmsh point: the sole FacetedVertex to reconstruct.
    node_group[n0.value] = wall;
    node_entity[n0.value] = 10;

    // e0 belongs to "Wall" (spanning vertex+curve+surface), e1..e3 belong to "Edges" as 3
    // separate 1-edge curves.
    auto e0 = mesh.add_edge(n0, n1);
    edge_group_var[e0.value] = wall;
    edge_entity_var[e0.value] = 20;
    auto e1 = mesh.add_edge(n1, n2);
    edge_group_var[e1.value] = edges_group;
    edge_entity_var[e1.value] = 21;
    auto e2 = mesh.add_edge(n2, n3);
    edge_group_var[e2.value] = edges_group;
    edge_entity_var[e2.value] = 22;
    auto e3 = mesh.add_edge(n3, n0);
    edge_group_var[e3.value] = edges_group;
    edge_entity_var[e3.value] = 23;

    // f0 belongs to "Wall" too, f1 to its own "Surfaces" group.
    auto f0 = mesh.add_face(n0, n1, n2);
    face_group_var[f0.value] = wall;
    face_entity_var[f0.value] = 30;
    auto f1 = mesh.add_face(n0, n2, n3);
    face_group_var[f1.value] = surfaces_group;
    face_entity_var[f1.value] = 31;

    const auto path = fixture_path("gecko_faceted_geometry_boundary_test.msh").string();
    io::SimplicialMeshWriter::write(path, mesh, groups);

    FacetedGeometry geom(path);
    std::filesystem::remove(path);

    SECTION("Entity counts: no volume for a triangles-only file") {
        REQUIRE(geom.nb_vertices() == 1);
        REQUIRE(geom.nb_curves() == 4);
        REQUIRE(geom.nb_surfaces() == 2);
        REQUIRE(geom.nb_volumes() == 0);
    }

    SECTION("A group spanning several dimensions returns entities of all of them") {
        const auto wall_entities = geom.entities(wall);
        REQUIRE(wall_entities.size() == 3);

        bool has_vertex = false, has_curve = false, has_surface = false;
        for (const auto &ref : wall_entities) {
            std::visit(
                [&](auto *entity) {
                    using T = std::decay_t<decltype(*entity)>;
                    if constexpr (std::is_same_v<T, FacetedVertex>) {
                        has_vertex = (entity->entity_tag() == 10);
                    } else if constexpr (std::is_same_v<T, FacetedCurve>) {
                        has_curve = (entity->entity_tag() == 20);
                    } else if constexpr (std::is_same_v<T, FacetedSurface>) {
                        has_surface = (entity->entity_tag() == 30);
                    }
                },
                ref);
        }
        REQUIRE(has_vertex);
        REQUIRE(has_curve);
        REQUIRE(has_surface);
    }

    SECTION("A single-dimension group returns only its own entities") {
        const auto surf_entities = geom.entities(surfaces_group);
        REQUIRE(surf_entities.size() == 1);
        REQUIRE(std::holds_alternative<const FacetedSurface *>(surf_entities[0]));
        REQUIRE(std::get<const FacetedSurface *>(surf_entities[0])->entity_tag() == 31);

        const auto edge_entities = geom.entities(edges_group);
        REQUIRE(edge_entities.size() == 3);
        for (const auto &ref : edge_entities) {
            REQUIRE(std::holds_alternative<const FacetedCurve *>(ref));
        }
    }

    SECTION("groups(dim) filters the model's physical groups by dimension") {
        auto dim1_groups = geom.groups(GroupDim::Dim1);
        std::size_t count = 0;
        for (const auto &g : dim1_groups) {
            REQUIRE(g.name == "Edges");
            ++count;
        }
        REQUIRE(count == 1);
    }

    SECTION("Distance, closest_point and project against a known FacetedSurface") {
        const auto &surface = find_by_tag(geom.surfaces(), 30); // triangle n0(0,0,0) n1(1,0,0) n2(1,1,0)
        const Point3d above_centroid(0.5, 0.25, 4.0);
        const Point3d cp = surface.closest_point(above_centroid);
        REQUIRE(cp.z() == Approx(0.0).margin(1e-9));
        REQUIRE(surface.distance(above_centroid) == Approx(4.0).margin(1e-9));

        Point3d mutated = above_centroid;
        surface.project(mutated);
        REQUIRE(mutated == cp);
    }

    SECTION("Entities can be looked up by their entity_tag") {
        const auto *vertex = geom.vertex_by_tag(10);
        REQUIRE(vertex != nullptr);
        REQUIRE(vertex->entity_tag() == 10);

        const auto *curve = geom.curve_by_tag(20);
        REQUIRE(curve != nullptr);
        REQUIRE(curve->entity_tag() == 20);

        const auto *surface = geom.surface_by_tag(30);
        REQUIRE(surface != nullptr);
        REQUIRE(surface->entity_tag() == 30);

        REQUIRE(geom.vertex_by_tag(999) == nullptr);
        REQUIRE(geom.curve_by_tag(999) == nullptr);
        REQUIRE(geom.surface_by_tag(999) == nullptr);
    }

    SECTION("Distance, closest_point and project against a known FacetedCurve") {
        const auto &curve = find_by_tag(geom.curves(), 20); // segment n0(0,0,0) -> n1(1,0,0)
        const Point3d p(0.5, 3.0, 0.0);
        const Point3d cp = curve.closest_point(p);
        REQUIRE(cp.x() == Approx(0.5).margin(1e-9));
        REQUIRE(cp.y() == Approx(0.0).margin(1e-9));
        REQUIRE(curve.distance(p) == Approx(3.0).margin(1e-9));

        Point3d mutated = p;
        curve.project(mutated);
        REQUIRE(mutated == cp);
    }
}

TEST_CASE("FacetedGeometry_AddsVolume_FromTetrahedra", "[FacetedGeometry]") {
    SimplicialMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(0, 1, 0);
    auto n3 = mesh.add_node(0, 0, 1);

    GroupRegistry groups;
    auto solid = groups.add_group("Solid", GroupDim::Dim3);

    auto &cell_group_var = mesh.add_variable<GroupId, CellType::Cell>(std::string(io::PHYSICAL_GROUP_VARIABLE));
    auto &cell_entity_var = mesh.add_variable<Int, CellType::Cell>(std::string(io::ENTITY_TAG_VARIABLE));

    auto c0 = mesh.add_cell(n0, n1, n2, n3);
    cell_group_var[c0.value] = solid;
    cell_entity_var[c0.value] = 40;

    const auto path = fixture_path("gecko_faceted_geometry_volume_test.msh").string();
    io::SimplicialMeshWriter::write(path, mesh, groups);

    FacetedGeometry geom(path);
    std::filesystem::remove(path);

    REQUIRE(geom.nb_volumes() == 1);
    REQUIRE(geom.volumes()[0].entity_tag() == 40);
    REQUIRE(geom.volumes()[0].dimension() == GroupDim::Dim3);
    REQUIRE(geom.volumes()[0].distance(Point3d(5.0, 5.0, 5.0)) == Approx(0.0).margin(1e-12));

    REQUIRE(geom.volume_by_tag(40) == &geom.volumes()[0]);
    REQUIRE(geom.volume_by_tag(999) == nullptr);

    const auto solid_entities = geom.entities(solid);
    REQUIRE(solid_entities.size() == 1);
    REQUIRE(std::holds_alternative<const FacetedVolume *>(solid_entities[0]));
}

TEST_CASE("FacetedGeometry_FromGmsh2DFile", "[FacetedGeometry]") {

    std::string dir(TEST_SAMPLES_DIR);
    const auto path = dir + "/square.msh";

    FacetedGeometry geom(path);

    REQUIRE(geom.nb_volumes() == 0);
    REQUIRE(geom.nb_surfaces() == 2);
    REQUIRE(geom.nb_curves() == 7);
    REQUIRE(geom.nb_vertices() == 6);

    //check the vertex tags
    auto tag_view_vertices =
        geom.vertices() | std::views::transform([](const FacetedVertex &v) { return v.entity_tag(); });
    REQUIRE_THAT(tag_view_vertices, Catch::Matchers::RangeEquals(std::vector{1, 2, 3, 4, 5, 6}));
    //check the curve tags
    auto tag_view_curves = geom.curves() | std::views::transform([](const FacetedCurve &c) { return c.entity_tag(); });
    REQUIRE_THAT(tag_view_curves, Catch::Matchers::RangeEquals(std::vector{1, 2, 3, 4, 5, 6, 7}));
    //check the surface tags
    auto tag_view_surfs =
        geom.surfaces() | std::views::transform([](const FacetedSurface &s) { return s.entity_tag(); });
    REQUIRE_THAT(tag_view_surfs, Catch::Matchers::RangeEquals(std::vector{1, 2}));

    REQUIRE(geom.surface_by_tag(1)->distance(Point3d(0.5, 0.5, 1.0)) == Approx(1.0).margin(1e-12));
    REQUIRE(geom.surface_by_tag(2)->distance(Point3d(1.5, 0.5, 1.0)) == Approx(1.0).margin(1e-12));
}

TEST_CASE("FacetedGeometry_FromGmsh3DFile", "[FacetedGeometry]") {

    std::string dir(TEST_SAMPLES_DIR);
    const auto path = dir + "/two_cubes.msh";

    FacetedGeometry geom(path);

    REQUIRE(geom.nb_volumes() == 2);
    REQUIRE(geom.nb_surfaces() == 11);
    REQUIRE(geom.nb_curves() == 20);
    REQUIRE(geom.nb_vertices() == 12);

    //check the vertex tags
    auto tag_view_vertices =
        geom.vertices() | std::views::transform([](const FacetedVertex &v) { return v.entity_tag(); });
    REQUIRE_THAT(tag_view_vertices, Catch::Matchers::RangeEquals(std::vector{1, 2, 3, 4, 5, 6, 7, 8, 13, 14, 15, 16}));
    //check the curve tags
    auto tag_view_curves = geom.curves() | std::views::transform([](const FacetedCurve &c) { return c.entity_tag(); });
    REQUIRE_THAT(tag_view_curves, Catch::Matchers::RangeEquals(std::vector{1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
                                                                           11, 12, 21, 22, 23, 24, 25, 26, 27, 28}));
    //check the surface tags
    auto tag_view_surfs =
        geom.surfaces() | std::views::transform([](const FacetedSurface &s) { return s.entity_tag(); });
    REQUIRE_THAT(tag_view_surfs, Catch::Matchers::RangeEquals(std::vector{1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 16}));

    //check the volume tags
    auto tag_view_volumes =
        geom.volumes() | std::views::transform([](const FacetedVolume &v) { return v.entity_tag(); });
    REQUIRE_THAT(tag_view_volumes, Catch::Matchers::RangeEquals(std::vector{1, 3}));
}

TEST_CASE("FacetedGeometry_FromGmshCylinderFile", "[FacetedGeometry]") {
    // test_data/cylinder.msh: an OpenCASCADE-generated tet mesh of a radius-0.5, length-2 cylinder
    // lying along x (Cylinder(1) = {0,0,0, 2,0,0, 0.5}), with no $PhysicalNames (like two_cubes.msh),
    // exercising entity-tag-only classification for a single-volume, curved-surface model: 2
    // vertices (the seam's 2 endpoints), 3 curves (the seam line + the 2 circular rims), 3 surfaces
    // (the lateral cylindrical surface + the 2 flat end caps) and 1 volume.
    std::string dir(TEST_SAMPLES_DIR);
    const auto path = dir + "/cylinder.msh";

    FacetedGeometry geom(path);

    REQUIRE(geom.nb_vertices() == 2);
    REQUIRE(geom.nb_curves() == 3);
    REQUIRE(geom.nb_surfaces() == 3);
    REQUIRE(geom.nb_volumes() == 1);

    auto tag_view_vertices =
        geom.vertices() | std::views::transform([](const FacetedVertex &v) { return v.entity_tag(); });
    REQUIRE_THAT(tag_view_vertices, Catch::Matchers::RangeEquals(std::vector{1, 2}));
    auto tag_view_curves = geom.curves() | std::views::transform([](const FacetedCurve &c) { return c.entity_tag(); });
    REQUIRE_THAT(tag_view_curves, Catch::Matchers::RangeEquals(std::vector{1, 2, 3}));
    auto tag_view_surfs =
        geom.surfaces() | std::views::transform([](const FacetedSurface &s) { return s.entity_tag(); });
    REQUIRE_THAT(tag_view_surfs, Catch::Matchers::RangeEquals(std::vector{1, 2, 3}));
    auto tag_view_volumes =
        geom.volumes() | std::views::transform([](const FacetedVolume &v) { return v.entity_tag(); });
    REQUIRE_THAT(tag_view_volumes, Catch::Matchers::RangeEquals(std::vector{1}));

    // The seam's 2 endpoints sit on the rim circles, at (x=2,y=0,z=0.5) and (x=0,y=0,z=0.5) — i.e.
    // radius 0.5 from the x axis at each end of the cylinder's length.
    const auto *v1 = geom.vertex_by_tag(1);
    REQUIRE(v1 != nullptr);
    REQUIRE(v1->closest_point(Point3d(0, 0, 0)).x() == Approx(2.0).margin(1e-9));
    REQUIRE(v1->closest_point(Point3d(0, 0, 0)).y() == Approx(0.0).margin(1e-9));
    REQUIRE(v1->closest_point(Point3d(0, 0, 0)).z() == Approx(0.5).margin(1e-9));
    const auto *v2 = geom.vertex_by_tag(2);
    REQUIRE(v2 != nullptr);
    REQUIRE(v2->closest_point(Point3d(0, 0, 0)).x() == Approx(0.0).margin(1e-9));

    // Surface 1 (the lateral surface) must be markedly closer to a point on the cylinder's own axis,
    // mid-length, than the flat end-cap surfaces 2 (at x=2) / 3 (at x=0) are.
    const auto *lateral = geom.surface_by_tag(1);
    const auto *cap_far = geom.surface_by_tag(2);
    const auto *cap_near = geom.surface_by_tag(3);
    REQUIRE(lateral != nullptr);
    REQUIRE(cap_far != nullptr);
    REQUIRE(cap_near != nullptr);
    const Point3d axis_midpoint(1.0, 0.0, 0.0);
    // The lateral surface's own triangulated facets slightly undershoot the analytic radius (a
    // chord-vs-arc faceting error), hence the looser margin than the flat, exactly-planar caps.
    REQUIRE(lateral->distance(axis_midpoint) == Approx(0.5).margin(5e-3));
    REQUIRE(cap_near->distance(axis_midpoint) == Approx(1.0).margin(1e-6));
    REQUIRE(cap_far->distance(axis_midpoint) == Approx(1.0).margin(1e-6));

    // The single volume's stub distance()/closest_point() (see FacetedVolume's own doc comment) is
    // always 0/identity, regardless of query point — still worth asserting explicitly here, since a
    // multi-surface curved model is exactly the case that stub is a known simplification for.
    const auto *volume = geom.volume_by_tag(1);
    REQUIRE(volume != nullptr);
    REQUIRE(volume->distance(axis_midpoint) == Approx(0.0).margin(1e-12));
}

TEST_CASE("FacetedSurface_Normal_ComesFromTheNearestFacet", "[FacetedGeometry]") {
    // 2 triangles meeting along a ridge, at right angles: one lies in z=0, the other in y=0. The
    // normal a query gets is the one of whichever facet it is nearest, which is all a faceted
    // surface has to give.
    SimplicialMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(1, 1, 0);
    auto n3 = mesh.add_node(1, 0, 1);

    std::vector<FaceId> faces{mesh.add_face(n0, n1, n2), mesh.add_face(n0, n1, n3)};
    const FacetedSurface surface(&mesh, faces, 1);

    // Just above the flat triangle, and nearer it than the upright one: its normal is along z.
    // Far above would not do — the upright triangle rises to z = 1 and reaches up to meet it.
    const Vector3d flat = surface.normal(Point3d(0.7, 0.3, 0.05));
    REQUIRE(std::abs(flat.z()) == Approx(1.0).margin(1e-12));

    // And just off the upright one: its normal is along y.
    const Vector3d upright = surface.normal(Point3d(0.9, -0.05, 0.5));
    REQUIRE(std::abs(upright.y()) == Approx(1.0).margin(1e-12));
}

TEST_CASE("FacetedSurface_PlaneSection_IsTheCurveThePlaneCutsOut", "[FacetedGeometry]") {
    // One square of 2 triangles in z=0, spanning [0,1]^2. The plane y = 1/2 cuts it along the
    // segment from (0,1/2,0) to (1,1/2,0), and that segment is what the section query answers on.
    SimplicialMesh mesh;
    auto n0 = mesh.add_node(0, 0, 0);
    auto n1 = mesh.add_node(1, 0, 0);
    auto n2 = mesh.add_node(1, 1, 0);
    auto n3 = mesh.add_node(0, 1, 0);
    std::vector<FaceId> faces{mesh.add_face(n0, n1, n2), mesh.add_face(n0, n2, n3)};
    const FacetedSurface surface(&mesh, faces, 1);

    const Point3d origin(0.0, 0.5, 0.0);
    const Vector3d normal(0.0, 1.0, 0.0);

    // A point above the middle of the section lands on it.
    const auto middle = surface.closest_point_on_section(Point3d(0.5, 0.5, 3.0), origin, normal);
    REQUIRE(middle.has_value());
    REQUIRE(middle->x() == Approx(0.5));
    REQUIRE(middle->y() == Approx(0.5));
    REQUIRE(middle->z() == Approx(0.0));

    // A point off to one side lands on the section's end, not on the nearest point of the surface:
    // the section stops where the surface does.
    const auto beyond = surface.closest_point_on_section(Point3d(4.0, 0.2, 0.0), origin, normal);
    REQUIRE(beyond.has_value());
    REQUIRE(beyond->x() == Approx(1.0));
    REQUIRE(beyond->y() == Approx(0.5));

    // A plane that misses the surface has no section at all, and says so rather than answering
    // with the nearest point of the surface — a caller has to be able to tell the difference.
    REQUIRE_FALSE(surface.closest_point_on_section(Point3d(0.5, 0.5, 0.0), Point3d(0.0, 5.0, 0.0), normal).has_value());
}
