#include <algorithm>
#include <filesystem>
#include <ranges>
#include <type_traits>
#include <variant>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <gecko/geom/FacetedGeometry.h>
#include <gecko/io/GmshMeshWriter.h>

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

    const auto solid_entities = geom.entities(solid);
    REQUIRE(solid_entities.size() == 1);
    REQUIRE(std::holds_alternative<const FacetedVolume *>(solid_entities[0]));
}
