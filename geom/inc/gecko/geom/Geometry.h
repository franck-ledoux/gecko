#pragma once

#include <vector>
#include <memory>
#include <functional>
#include <cstdint>

#include <gecko/geom/GeomEntities.h>
#include <gecko/geom/GeomEntityRefVariant.h>
#include <gecko/utils/Types.h>
#include <gecko/math/Point3d.h>

namespace gecko {

    // Identifiants forts distincts pour chaque entité CAO
    struct GeomVertexTag {};
    using GeomVertexId = StrongId<GeomVertexTag, UInt>;
    struct GeomCurveTag {};
    using GeomCurveId = StrongId<GeomCurveTag, UInt>;
    struct GeomSurfaceTag {};
    using GeomSurfaceId = StrongId<GeomSurfaceTag, UInt>;
    struct GeomVolumeTag {};
    using GeomVolumeId = StrongId<GeomVolumeTag, UInt>;

    /**
     * @class Geometry
     * @brief Owns the CAD model, storing all its vertices, curves, surfaces and volumes,
     * and providing strongly-typed id-based access to them.
     */
    class Geometry {
    public:
        /** @brief Default constructor. Builds an empty geometric model. */
        Geometry() = default;
        ~Geometry() = default;

        // Désactivation de la copie pour éviter les réallocations non contrôlées de la géométrie
        /** @brief Copy is disabled to avoid uncontrolled reallocations of the geometry. */
        Geometry(const Geometry &) = delete;
        /** @brief Copy is disabled to avoid uncontrolled reallocations of the geometry. @return Reference to this Geometry. */
        Geometry &operator=(const Geometry &) = delete;

        // Autorisation du déplacement (Move semantics)
        /** @brief Move constructor. */
        Geometry(Geometry &&) noexcept = default;
        /** @brief Move assignment. @return Reference to this Geometry. */
        Geometry &operator=(Geometry &&) noexcept = default;

        // ===================================================================
        // 1. AJOUT D'ENTITÉS
        // ===================================================================

        /**
         * @brief Adds a new vertex to the model.
         * @param position Position of the vertex.
         * @return The id of the newly created vertex.
         */
        GeomVertexId add_vertex(Point3d position) {
            auto id = GeomVertexId{static_cast<std::uint32_t>(m_vertices.size())};
            //m_vertices.emplace_back(position);
            return id;
        }

        /**
         * @brief Adds a new curve to the model.
         * @param length Length parameter of the curve (default: 1.0).
         * @return The id of the newly created curve.
         */
        GeomCurveId add_curve(double length = 1.0) {
            auto id = GeomCurveId{static_cast<std::uint32_t>(m_curves.size())};
            //m_curves.emplace_back(length);
            return id;
        }

        /**
         * @brief Adds a new surface to the model.
         * @param area Area parameter of the surface (default: 1.0).
         * @return The id of the newly created surface.
         */
        GeomSurfaceId add_surface(double area = 1.0) {
            auto id = GeomSurfaceId{static_cast<std::uint32_t>(m_surfaces.size())};
            //m_surfaces.emplace_back(area);
            return id;
        }

        /**
         * @brief Adds a new volume to the model.
         * @param volume Volume parameter (default: 1.0).
         * @return The id of the newly created volume.
         */
        GeomVolumeId add_volume(double volume = 1.0) {
            auto id = GeomVolumeId{static_cast<std::uint32_t>(m_volumes.size())};
            //m_volumes.emplace_back(volume);
            return id;
        }

        // ===================================================================
        // 2. ACCÈS DIRECT AUX ENTITÉS PAR ID (Accès typé fort)
        // ===================================================================

        /**
         * @brief Accesses a vertex by its id.
         * @param id Id of the vertex to retrieve.
         * @return Const reference to the requested vertex.
         * @pre @p id must be a valid id previously returned by add_vertex().
         */
        [[nodiscard]] const GeomVertex &get_vertex(GeomVertexId id) const { return m_vertices[id.value]; }

        /**
         * @brief Accesses a curve by its id.
         * @param id Id of the curve to retrieve.
         * @return Const reference to the requested curve.
         * @pre @p id must be a valid id previously returned by add_curve().
         */
        [[nodiscard]] const GeomCurve &get_curve(GeomCurveId id) const { return m_curves[id.value]; }

        /**
         * @brief Accesses a surface by its id.
         * @param id Id of the surface to retrieve.
         * @return Const reference to the requested surface.
         * @pre @p id must be a valid id previously returned by add_surface().
         */
        [[nodiscard]] const GeomSurface &get_surface(GeomSurfaceId id) const { return m_surfaces[id.value]; }

        /**
         * @brief Accesses a volume by its id.
         * @param id Id of the volume to retrieve.
         * @return Const reference to the requested volume.
         * @pre @p id must be a valid id previously returned by add_volume().
         */
        [[nodiscard]] const GeomVolume &get_volume(GeomVolumeId id) const { return m_volumes[id.value]; }

        // ===================================================================
        // 3. OBTENTION DE LA RÉFÉRENCE VARIANT POUR LE MAILLAGE / LISSEUR
        // ===================================================================

        /**
         * @brief Gets a type-erased reference to a vertex, usable uniformly with other entity kinds.
         * @param id Id of the vertex.
         * @return A GeomEntityRefVariant wrapping a pointer to the vertex.
         */
        [[nodiscard]] GeomEntityRefVariant get_entity_ref(GeomVertexId id) const { return &get_vertex(id); }

        /**
         * @brief Gets a type-erased reference to a curve, usable uniformly with other entity kinds.
         * @param id Id of the curve.
         * @return A GeomEntityRefVariant wrapping a pointer to the curve.
         */
        [[nodiscard]] GeomEntityRefVariant get_entity_ref(GeomCurveId id) const { return &get_curve(id); }

        /**
         * @brief Gets a type-erased reference to a surface, usable uniformly with other entity kinds.
         * @param id Id of the surface.
         * @return A GeomEntityRefVariant wrapping a pointer to the surface.
         */
        [[nodiscard]] GeomEntityRefVariant get_entity_ref(GeomSurfaceId id) const { return &get_surface(id); }

        /**
         * @brief Gets a type-erased reference to a volume, usable uniformly with other entity kinds.
         * @param id Id of the volume.
         * @return A GeomEntityRefVariant wrapping a pointer to the volume.
         */
        [[nodiscard]] GeomEntityRefVariant get_entity_ref(GeomVolumeId id) const { return &get_volume(id); }

        // ===================================================================
        // 4. COMPTEURS ET VÉRIFICATIONS
        // ===================================================================

        /** @brief Gets the number of vertices in the model. @return The vertex count. */
        [[nodiscard]] std::size_t vertex_count() const noexcept { return m_vertices.size(); }
        /** @brief Gets the number of curves in the model. @return The curve count. */
        [[nodiscard]] std::size_t curve_count() const noexcept { return m_curves.size(); }
        /** @brief Gets the number of surfaces in the model. @return The surface count. */
        [[nodiscard]] std::size_t surface_count() const noexcept { return m_surfaces.size(); }
        /** @brief Gets the number of volumes in the model. @return The volume count. */
        [[nodiscard]] std::size_t volume_count() const noexcept { return m_volumes.size(); }

    private:
        std::vector<GeomVertex> m_vertices;
        std::vector<GeomCurve> m_curves;
        std::vector<GeomSurface> m_surfaces;
        std::vector<GeomVolume> m_volumes;
    };

} // namespace gecko