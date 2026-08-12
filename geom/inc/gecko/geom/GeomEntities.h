#pragma once

#include <gecko/math/Point3d.h>
#include <gecko/math/Vector3d.h>
#include <gecko/utils/Groups.h>
#include <gecko/geom_itf/Concepts.h>

namespace gecko {

    /**
     * @class GeomVertex
     * @brief CAD entity representing a fixed point (0D) in the geometric model.
     */
    class GeomVertex {
    public:
        /**
         * @brief Constructor.
         * @param pos Position of the vertex.
         */
        explicit GeomVertex(Point3d pos) : m_pos(pos) {}

        /**
         * @brief Projects any point onto this vertex.
         * @return The vertex's own position, regardless of the input point.
         */
        [[nodiscard]] Point3d project(const Point3d &) const noexcept { return m_pos; }

        /**
         * @brief Computes the distance between a point and this vertex.
         * @param p Query point.
         * @return The Euclidean distance between @p p and the vertex's position.
         */
        [[nodiscard]] double distance(const Point3d &p) const noexcept { return Vector3d(m_pos, p).norm(); }

        /** @brief Gets the topological dimension of this entity. @return GroupDim::Dim0. */
        [[nodiscard]] constexpr GroupDim dimension() const noexcept { return GroupDim::Dim0; }

    private:
        Point3d m_pos;
    };
    static_assert(GeomEntityConcept<GeomVertex>, "GeomVertex must satisfy GeomEntityConcept");

    // 2. Courbe CAO (1D)
    /**
     * @class GeomCurve
     * @brief CAD entity representing a curve (1D) in the geometric model.
     */
    class GeomCurve {
    public:
        /**
         * @brief Projects a point onto this curve.
         * @param p Query point.
         * @return The closest point of the curve to @p p.
         */
        [[nodiscard]] Point3d project(const Point3d &p) const {
            return p; // Exemple
        }

        /**
         * @brief Computes the distance between a point and this curve.
         * @param p Query point.
         * @return The Euclidean distance between @p p and its projection onto the curve.
         */
        [[nodiscard]] double distance(const Point3d &p) const { return Vector3d(project(p), p).norm(); }

        /** @brief Gets the topological dimension of this entity. @return GroupDim::Dim1. */
        [[nodiscard]] constexpr GroupDim dimension() const noexcept { return GroupDim::Dim1; }
    };
    static_assert(GeomEntityConcept<GeomCurve>, "GeomCurve must satisfy GeomEntityConcept");

    /**
     * @class GeomSurface
     * @brief CAD entity representing a surface (2D) in the geometric model.
     */
    class GeomSurface {
    public:
        /**
         * @brief Projects a point onto this surface.
         * @param p Query point.
         * @return The closest point of the surface to @p p.
         */
        [[nodiscard]] Point3d project(const Point3d &p) const {
            // Logique de projection 2D (ex: projection sur une NURBS)
            return p; // Exemple
        }

        /**
         * @brief Computes the distance between a point and this surface.
         * @param p Query point.
         * @return The Euclidean distance between @p p and its projection onto the surface.
         */
        [[nodiscard]] double distance(const Point3d &p) const { return Vector3d(project(p), p).norm(); }

        /** @brief Gets the topological dimension of this entity. @return GroupDim::Dim2. */
        [[nodiscard]] constexpr GroupDim dimension() const noexcept { return GroupDim::Dim2; }
    };
    static_assert(GeomEntityConcept<GeomSurface>, "GeomSurface must satisfy GeomEntityConcept");

    // 4. Volume CAO (3D)
    /**
     * @class GeomVolume
     * @brief CAD entity representing a volume (3D) in the geometric model.
     */
    class GeomVolume {
    public:
        /**
         * @brief Projects a point onto this volume.
         * @param p Query point.
         * @return The input point unchanged, since any point already lies within (or maps freely onto) the volume.
         */
        [[nodiscard]] Point3d project(const Point3d &p) const noexcept { return p; } // Déplacement libre

        /**
         * @brief Computes the distance between a point and this volume.
         * @return Always 0.0, since points are considered free to move within the volume.
         */
        [[nodiscard]] double distance(const Point3d & /*p*/) const noexcept { return 0.0; }

        /** @brief Gets the topological dimension of this entity. @return GroupDim::Dim3. */
        [[nodiscard]] constexpr GroupDim dimension() const noexcept { return GroupDim::Dim3; }
    };
    static_assert(GeomEntityConcept<GeomVolume>, "GeomVolume must satisfy GeomEntityConcept");

} // namespace gecko