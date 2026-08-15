#pragma once

#include <concepts>
#include <ranges>
#include <utility>
#include <gecko/math/Point3d.h>
#include <gecko/utils/Groups.h>

namespace gecko {

    /**
     * @brief Concept that must be satisfied by geom entities (vertices, curves, surfaces,
     * volumes) to be used for block and mesh operations.
     *
     * A conforming type must provide: a distance() method returning the distance to a query
     * point; a closest_point() method returning the closest point of the entity to a query point
     * (a pure query, does not modify its argument); a project() method that snaps a point onto
     * the entity **in place** (mutates its argument, returns nothing — unrelated to
     * BezierCurve::project(), which instead returns a curvilinear parameter; easy to conflate
     * given the shared name, but a different type in a different module); and a dimension()
     * method returning the entity's topological dimension.
     * @tparam T Candidate geometric entity type.
     */
    template<typename T>
    concept GeomEntityConcept = requires(const T entity, const Point3d &p, Point3d &mutable_p) {
        { entity.distance(p) } -> std::same_as<double>;
        { entity.closest_point(p) } -> std::same_as<Point3d>;
        { entity.project(mutable_p) } -> std::same_as<void>;
        { entity.dimension() } -> std::same_as<GroupDim>;
    };

    /**
     * @brief Concept describing a geometric model: a collection of GeomEntityConcept-conforming
     * vertices, curves, surfaces and volumes, plus the model's named groups of entities (e.g.
     * from Gmsh `$PhysicalNames`), which may span several dimensions at once (a "Wall" group
     * could contain 1 surface, 4 curves and 4 vertices together).
     *
     * A conforming type must provide: nb_vertices()/nb_curves()/nb_surfaces()/nb_volumes() entity
     * counts; vertices()/curves()/surfaces()/volumes() ranges of GeomEntityConcept-conforming
     * entities (any range type, not necessarily a `std::vector`); groups() (every group,
     * regardless of dimension) and groups(GroupDim) (groups of one dimension only); and
     * entities(GroupId), the (possibly dimension-mixed) entities belonging to a given group.
     * @tparam T Candidate geometric model type.
     */
    template<typename T>
    concept GeomModelConcept =
        requires(const T model, GroupId gid, GroupDim dim) {
            { model.nb_vertices() } -> std::convertible_to<std::size_t>;
            { model.nb_curves() } -> std::convertible_to<std::size_t>;
            { model.nb_surfaces() } -> std::convertible_to<std::size_t>;
            { model.nb_volumes() } -> std::convertible_to<std::size_t>;
            { model.vertices() } -> std::ranges::range;
            { model.curves() } -> std::ranges::range;
            { model.surfaces() } -> std::ranges::range;
            { model.volumes() } -> std::ranges::range;
            { model.groups() } -> std::ranges::range;
            { model.groups(dim) } -> std::ranges::range;
            { model.entities(gid) } -> std::ranges::range;
        } && GeomEntityConcept<std::ranges::range_value_t<decltype(std::declval<const T &>().vertices())>> &&
        GeomEntityConcept<std::ranges::range_value_t<decltype(std::declval<const T &>().curves())>> &&
        GeomEntityConcept<std::ranges::range_value_t<decltype(std::declval<const T &>().surfaces())>> &&
        GeomEntityConcept<std::ranges::range_value_t<decltype(std::declval<const T &>().volumes())>>;

} // namespace gecko