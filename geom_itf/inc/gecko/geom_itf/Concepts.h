#pragma once

#include <concepts>
#include <gecko/math/Point3d.h>
#include <gecko/utils/Groups.h>

namespace gecko {

    /**
     * @brief Concept that must be satisfied by geom entities to be used for block and mesh operations.
     *
     * A conforming type must provide a project() method returning the closest Point3d,
     * a distance() method returning the distance to a query point, and a dimension()
     * method returning its topological dimension.
     * @tparam T Candidate geometric entity type.
     */
    template<typename T>
    concept GeomEntityConcept = requires(const T entity, const Point3d &p) {
        { entity.project(p) } -> std::same_as<Point3d>;
        { entity.distance(p) } -> std::same_as<double>;
        { entity.dimension() } -> std::same_as<GroupDim>;
    };

    /**
     * @brief Concept describing any geometric surface (faceted, analytical, CAD, etc.).
     *
     * A conforming type must provide project(), normal_at(), distance_to() and tag() methods.
     * @tparam T Candidate surface type.
     */
    // Concept décrivant n'importe quelle surface géométrique (facettée, analytique, CAD, etc.)
    template<typename T>
    concept GeometryFace = requires(const T face, const Point3d &p) {
        // Exige une méthode project() retournant un Point3d
        { face.project(p) } -> std::same_as<Point3d>;

        // Exige une méthode normal_at() retournant un Vector3d
        { face.normal_at(p) } -> std::same_as<Vector3d>;

        // Exige une méthode distance_to() retournant un double
        { face.distance_to(p) } -> std::convertible_to<double>;

        // Exige un tag/identifiant
        { face.tag() } -> std::integral;
    };

    /**
     * @brief Concept describing a 1D geometric curve.
     *
     * A conforming type must provide project(), point_at() (parameterized in [0, 1])
     * and tag() methods.
     * @tparam T Candidate curve type.
     */
    // Concept décrivant une courbe géométrique 1D
    template<typename T>
    concept GeometryEdge = requires(const T edge, const Point3d &p, double t) {
        { edge.project(p) } -> std::same_as<Point3d>;
        { edge.point_at(t) } -> std::same_as<Point3d>; // Paramétrage [0, 1]
        { edge.tag() } -> std::integral;
    };

} // namespace gecko