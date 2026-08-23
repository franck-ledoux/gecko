#pragma once

#include <concepts>
#include <cstddef>
#include <ranges>
#include <utility>
#include <gecko/math/BezierCurve.h>
#include <gecko/math/CurveSurfaceTraits.h>
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
     * given the shared name, but a different type in a different module); a dimension()
     * method returning the entity's topological dimension; and a tag that identify the cell in its dimension, using
     * Gmsh approach.
     * @tparam T Candidate geometric entity type.
     */
    template<typename T>
    concept GeomEntityConcept = requires(const T entity, const Point3d &p, Point3d &mutable_p) {
        { entity.distance(p) } -> std::same_as<double>;
        { entity.closest_point(p) } -> std::same_as<Point3d>;
        { entity.project(mutable_p) } -> std::same_as<void>;
        { entity.dimension() } -> std::same_as<GroupDim>;
        { entity.entity_tag() } -> std::same_as<Int>;
    };

    /**
     * @brief Concept describing a geometric model: a collection of GeomEntityConcept-conforming
     * vertices, curves, surfaces and volumes, plus the model's named groups of entities (e.g.
     * from Gmsh `$PhysicalNames`), which may span several dimensions at once (a "Wall" group
     * could contain 1 surface, 4 curves and 4 vertices together).
     *
     * A conforming type must provide: nb_vertices()/nb_curves()/nb_surfaces()/nb_volumes() entity
     * counts; vertices()/curves()/surfaces()/volumes() ranges of GeomEntityConcept-conforming
     * entities (any range type, not necessarily a `std::vector`); vertex_by_tag()/curve_by_tag()/
     * surface_by_tag()/volume_by_tag(), each looking up an entity by its entity_tag() (see
     * GeomEntityConcept) and returning a pointer to it, or `nullptr` if no entity of that
     * dimension has that tag; groups() (every group, regardless of dimension) and
     * groups(GroupDim) (groups of one dimension only); entities(GroupId), the (possibly
     * dimension-mixed) entities belonging to a given group; and containing_entities(GroupDim, Int),
     * the model's B-Rep incidence — an entity plus every higher-dimensional entity containing it,
     * as (dimension, entity_tag) pairs — which lets a block cell's classification be inferred from
     * that of its own boundary rather than from proximity alone (see `Blocking::classify()`).
     * @tparam T Candidate geometric model type.
     */
    template<typename T>
    concept GeomModelConcept =
        requires(const T model, GroupId gid, GroupDim dim, Int tag) {
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
            { model.containing_entities(dim, tag) } -> std::ranges::range;
            {
                model.vertex_by_tag(tag)
            } -> std::same_as<const std::ranges::range_value_t<decltype(std::declval<const T &>().vertices())> *>;
            {
                model.curve_by_tag(tag)
            } -> std::same_as<const std::ranges::range_value_t<decltype(std::declval<const T &>().curves())> *>;
            // Curves alone must also report a direction: fitting a curved block edge onto one needs
            // its tangent at the edge's ends, and unlike a point projection there is no meaningful
            // counterpart of this on a vertex, surface or volume.
            { model.curve_by_tag(tag)->tangent(std::declval<const Point3d &>()) } -> std::same_as<Vector3d>;
            {
                model.surface_by_tag(tag)
            } -> std::same_as<const std::ranges::range_value_t<decltype(std::declval<const T &>().surfaces())> *>;
            {
                model.volume_by_tag(tag)
            } -> std::same_as<const std::ranges::range_value_t<decltype(std::declval<const T &>().volumes())> *>;
        } && GeomEntityConcept<std::ranges::range_value_t<decltype(std::declval<const T &>().vertices())>> &&
        GeomEntityConcept<std::ranges::range_value_t<decltype(std::declval<const T &>().curves())>> &&
        GeomEntityConcept<std::ranges::range_value_t<decltype(std::declval<const T &>().surfaces())>> &&
        GeomEntityConcept<std::ranges::range_value_t<decltype(std::declval<const T &>().volumes())>>;

    /**
     * @brief Concept describing a control-point-based curve representation usable as a `block`
     * edge's geometry (e.g. `gecko::BezierCurve<N, Point3d>` today; a future NURBS/B-spline curve
     * satisfying the same concept plugs in without touching any code written against it).
     *
     * A conforming type must be default-constructible (needed to serve as a CGAL `Cell_attribute`
     * info payload) and provide: `Degree`, the curve's polynomial/parametric degree; `NumControlPoints`
     * (`Degree + 1`), the number of degrees of freedom driving the curve; `value(t)`, evaluating the
     * curve at a parametric coordinate in `[0.0, 1.0]`; `control_points()`, a zero-copy, ordered,
     * endpoint-inclusive view of its control points; indexed read/write access (`operator[]`) to
     * individual control points — used together to linearly initialize a fresh curve between two
     * corners, to reposition interior control points once classified onto a geometric curve, and to
     * build a reversed copy of a curve (needed when orienting a block's edges consistently around a
     * face/volume for Coons/TFI construction); and a `CurveSurfaceTraits<T>::Surface` specialization
     * (see `math/CurveSurfaceTraits.h`) naming the paired tensor-product surface representation
     * `Blocking`/`coons_surface_from_edges()` build faces with — the seam that keeps `block` generic
     * over the representation family instead of hardcoding `BezierSurface`/`BezierVolume`.
     * @tparam T Candidate curve representation type.
     */
    template<typename T>
    concept EdgeCurveConcept =
        std::default_initializable<T> && requires(T curve, const T const_curve, double t, std::size_t i) {
            // Degree is asked of the *object*, not of the type: a block structure has to be able to
            // change order while it is being edited, which a degree fixed in the type cannot do.
            { const_curve.degree() } -> std::convertible_to<std::size_t>;
            { const_curve.nb_control_points() } -> std::convertible_to<std::size_t>;
            { T(i) };
            { const_curve.value(t) } -> std::same_as<Point3d>;
            { const_curve.control_points() } -> std::ranges::range;
            { const_curve[i] } -> std::same_as<const Point3d &>;
            { curve[i] } -> std::same_as<Point3d &>;
            typename CurveSurfaceTraits<T>::Surface;
        };
    static_assert(EdgeCurveConcept<BezierCurve<Point3d>>, "BezierCurve<Point3d> must satisfy EdgeCurveConcept");

} // namespace gecko