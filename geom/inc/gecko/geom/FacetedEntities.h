#pragma once

#include <cassert>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include <gecko/math/ClosestPoint.h>
#include <gecko/math/Point3d.h>
#include <gecko/math/Vector3d.h>
#include <gecko/mesh/UnstructuredMesh.h>
#include <gecko/utils/Groups.h>
#include <gecko/utils/Types.h>
#include <gecko/geom_itf/Concepts.h>

namespace gecko {

    /**
     * @class FacetedVertex
     * @brief Faceted (mesh-backed) CAD vertex: a single node of a SimplicialMesh, reconstructed
     * from a Gmsh point element (type 15).
     *
     * Stores a non-owning pointer to the backing mesh plus the id of its node, rather than a
     * copied position, so it stays in sync with the mesh and avoids duplicating coordinates.
     */
    class FacetedVertex {
    public:
        /**
         * @brief Constructor.
         * @param mesh Backing mesh; must outlive this vertex.
         * @param node Id of the mesh node this vertex is built from.
         * @param entity_tag Gmsh elementary entity tag this vertex was reconstructed from.
         */
        FacetedVertex(const SimplicialMesh *mesh, NodeId node, Int entity_tag)
            : m_mesh(mesh), m_node(node), m_entity_tag(entity_tag) {}

        /**
         * @brief Gets the closest point of this vertex to a query point.
         * @return The vertex's own position, regardless of the input point.
         */
        [[nodiscard]] Point3d closest_point(const Point3d &) const { return m_mesh->node(m_node); }

        /**
         * @brief Snaps a point onto this vertex, in place.
         * @param p Point to snap; overwritten with the vertex's own position.
         */
        void project(Point3d &p) const { p = closest_point(p); }

        /**
         * @brief Computes the distance between a point and this vertex.
         * @param p Query point.
         * @return The Euclidean distance between @p p and the vertex's position.
         */
        [[nodiscard]] double distance(const Point3d &p) const { return Vector3d(m_mesh->node(m_node), p).norm(); }

        /** @brief Gets the topological dimension of this entity. @return GroupDim::Dim0. */
        [[nodiscard]] constexpr GroupDim dimension() const noexcept { return GroupDim::Dim0; }

        /** @brief Gets the Gmsh elementary entity tag this vertex was reconstructed from. @return The entity tag. */
        [[nodiscard]] Int entity_tag() const noexcept { return m_entity_tag; }

        /**
         * @brief Gets the backing mesh node this vertex was reconstructed from — the counterpart of
         * FacetedCurve::edges()/FacetedSurface::faces(), letting entity incidence be derived from
         * shared mesh cells (see FacetedGeometry::containing_entities()).
         * @return The mesh node id.
         */
        [[nodiscard]] NodeId node() const noexcept { return m_node; }

    private:
        const SimplicialMesh *m_mesh;
        NodeId m_node;
        Int m_entity_tag;
    };
    static_assert(GeomEntityConcept<FacetedVertex>, "FacetedVertex must satisfy GeomEntityConcept");

    /**
     * @class FacetedCurve
     * @brief Faceted (mesh-backed) CAD curve: the set of mesh edges sharing one Gmsh elementary
     * entity tag.
     *
     * Stores a non-owning pointer to the backing mesh plus the ids of its constituent edges,
     * rather than copied coordinates. Queries are brute-force over every edge (no acceleration
     * structure), the closest point of the segment they define to the query point.
     */
    class FacetedCurve {
    public:
        /**
         * @brief Constructor.
         * @param mesh Backing mesh; must outlive this curve.
         * @param edges Ids of the mesh edges this curve is built from (must be non-empty).
         * @param entity_tag Gmsh elementary entity tag this curve was reconstructed from.
         */
        FacetedCurve(const SimplicialMesh *mesh, std::vector<EdgeId> edges, Int entity_tag)
            : m_mesh(mesh), m_edges(std::move(edges)), m_entity_tag(entity_tag) {
            assert(!m_edges.empty() && "FacetedCurve: must be built from at least one edge");
        }

        /**
         * @brief Gets the closest point of this curve to a query point.
         * @param p Query point.
         * @return The closest point, over every constituent edge, of that edge's segment to @p p.
         */
        [[nodiscard]] Point3d closest_point(const Point3d &p) const {
            Point3d best{};
            double best_dist_sq = std::numeric_limits<double>::infinity();
            for (EdgeId e : m_edges) {
                const auto nodes = m_mesh->edge_nodes(e);
                const Point3d cand = closest_point_on_segment(p, m_mesh->node(nodes[0]), m_mesh->node(nodes[1]));
                if (const double d = Vector3d(cand, p).norm_sq(); d < best_dist_sq) {
                    best_dist_sq = d;
                    best = cand;
                }
            }
            return best;
        }

        /**
         * @brief Snaps a point onto this curve, in place.
         * @param p Point to snap; overwritten with its closest point on the curve.
         */
        void project(Point3d &p) const { p = closest_point(p); }

        /**
         * @brief Computes the distance between a point and this curve.
         * @param p Query point.
         * @return The Euclidean distance between @p p and its closest point on the curve.
         */
        [[nodiscard]] double distance(const Point3d &p) const { return Vector3d(closest_point(p), p).norm(); }

        /** @brief Gets the topological dimension of this entity. @return GroupDim::Dim1. */
        [[nodiscard]] constexpr GroupDim dimension() const noexcept { return GroupDim::Dim1; }

        /** @brief Gets the Gmsh elementary entity tag this curve was reconstructed from. @return The entity tag. */
        [[nodiscard]] Int entity_tag() const noexcept { return m_entity_tag; }

        /** @brief Gets a zero-copy view onto the ids of this curve's constituent mesh edges. @return A span of edge ids. */
        [[nodiscard]] std::span<const EdgeId> edges() const noexcept { return m_edges; }

    private:
        const SimplicialMesh *m_mesh;
        std::vector<EdgeId> m_edges;
        Int m_entity_tag;
    };
    static_assert(GeomEntityConcept<FacetedCurve>, "FacetedCurve must satisfy GeomEntityConcept");

    /**
     * @class FacetedSurface
     * @brief Faceted (mesh-backed) CAD surface: the set of mesh faces (triangles) sharing one
     * Gmsh elementary entity tag.
     *
     * Stores a non-owning pointer to the backing mesh plus the ids of its constituent faces,
     * rather than copied coordinates. Closest-point/distance queries are a brute-force scan over
     * those triangles via math::closest_point_on_triangle() — no spatial acceleration structure
     * (e.g. an AABB tree) is used, deliberately: real-world faceted surfaces in this project's
     * scope are small enough that a linear scan is fast enough, and it avoids depending on
     * CGAL's own (vendored, unpatched) template-heavy acceleration-structure headers, whose
     * exact grammar acceptance varies across compiler versions.
     */
    class FacetedSurface {
    public:
        /**
         * @brief Constructor. Caches @p faces' corner points for repeated closest-point queries.
         * @param mesh Backing mesh; must outlive this surface.
         * @param faces Ids of the mesh faces this surface is built from (must be non-empty).
         * @param entity_tag Gmsh elementary entity tag this surface was reconstructed from.
         */
        FacetedSurface(const SimplicialMesh *mesh, std::vector<FaceId> faces, Int entity_tag)
            : m_mesh(mesh),
              m_faces(std::move(faces)),
              m_entity_tag(entity_tag),
              m_triangles(build_triangles(m_mesh, m_faces)) {
            assert(!m_faces.empty() && "FacetedSurface: must be built from at least one face");
        }

        /**
         * @brief Gets the closest point of this surface to a query point.
         * @param p Query point.
         * @return The closest point of the surface (over every constituent triangle) to @p p.
         */
        [[nodiscard]] Point3d closest_point(const Point3d &p) const {
            Point3d best;
            double best_dist_sq = std::numeric_limits<double>::max();
            for (const auto &[a, b, c] : m_triangles) {
                const Point3d candidate = closest_point_on_triangle(p, a, b, c);
                if (const double dist_sq = Vector3d(p, candidate).norm_sq(); dist_sq < best_dist_sq) {
                    best_dist_sq = dist_sq;
                    best = candidate;
                }
            }
            return best;
        }

        /**
         * @brief Snaps a point onto this surface, in place.
         * @param p Point to snap; overwritten with its closest point on the surface.
         */
        void project(Point3d &p) const { p = closest_point(p); }

        /**
         * @brief Computes the distance between a point and this surface.
         * @param p Query point.
         * @return The Euclidean distance between @p p and its closest point on the surface.
         */
        [[nodiscard]] double distance(const Point3d &p) const { return Vector3d(closest_point(p), p).norm(); }

        /** @brief Gets the topological dimension of this entity. @return GroupDim::Dim2. */
        [[nodiscard]] constexpr GroupDim dimension() const noexcept { return GroupDim::Dim2; }

        /** @brief Gets the Gmsh elementary entity tag this surface was reconstructed from. @return The entity tag. */
        [[nodiscard]] Int entity_tag() const noexcept { return m_entity_tag; }

        /** @brief Gets a zero-copy view onto the ids of this surface's constituent mesh faces. @return A span of face ids. */
        [[nodiscard]] std::span<const FaceId> faces() const noexcept { return m_faces; }

    private:
        /** @brief One triangle's 3 corner points, cached so closest_point() doesn't re-hit the mesh. */
        struct Triangle {
            Point3d a, b, c;
        };

        /**
         * @brief Builds the list of triangles making up @p faces, for closest-point queries.
         * @param mesh Backing mesh to read node positions from.
         * @param faces Ids of the mesh faces to convert.
         * @return One Triangle per face, in the same order.
         */
        static std::vector<Triangle> build_triangles(const SimplicialMesh *mesh, std::span<const FaceId> faces) {
            std::vector<Triangle> triangles;
            triangles.reserve(faces.size());
            for (FaceId f : faces) {
                const auto nodes = mesh->face_nodes(f);
                triangles.push_back({mesh->node(nodes[0]), mesh->node(nodes[1]), mesh->node(nodes[2])});
            }
            return triangles;
        }

        const SimplicialMesh *m_mesh;
        std::vector<FaceId> m_faces;
        Int m_entity_tag;
        std::vector<Triangle> m_triangles;
    };
    static_assert(GeomEntityConcept<FacetedSurface>, "FacetedSurface must satisfy GeomEntityConcept");

    /**
     * @class FacetedVolume
     * @brief Faceted CAD volume reconstructed from a Gmsh elementary entity tag shared by tets.
     *
     * Real point-in-tet-mesh closest-point/distance queries are a separate, substantial problem
     * (out of scope for this pass); this class mirrors the existing GeomVolume stub's trivial
     * semantics instead (any point is considered to already lie within/on the volume).
     */
    class FacetedVolume {
    public:
        /**
         * @brief Constructor.
         * @param entity_tag Gmsh elementary entity tag this volume was reconstructed from.
         */
        explicit FacetedVolume(Int entity_tag) : m_entity_tag(entity_tag) {}

        /**
         * @brief Gets the closest point of this volume to a query point.
         * @param p Query point.
         * @return The input point unchanged, since any point already lies within (or maps freely onto) the volume.
         */
        [[nodiscard]] Point3d closest_point(const Point3d &p) const noexcept { return p; }

        /**
         * @brief Snaps a point onto this volume, in place.
         * @param p Point to snap; left unchanged, since any point already lies within (or maps freely onto) the volume.
         */
        void project(Point3d &p) const noexcept { p = closest_point(p); }

        /**
         * @brief Computes the distance between a point and this volume.
         * @return Always 0.0, since points are considered free to move within the volume.
         */
        [[nodiscard]] double distance(const Point3d & /*p*/) const noexcept { return 0.0; }

        /** @brief Gets the topological dimension of this entity. @return GroupDim::Dim3. */
        [[nodiscard]] constexpr GroupDim dimension() const noexcept { return GroupDim::Dim3; }

        /** @brief Gets the Gmsh elementary entity tag this volume was reconstructed from. @return The entity tag. */
        [[nodiscard]] Int entity_tag() const noexcept { return m_entity_tag; }

    private:
        Int m_entity_tag;
    };
    static_assert(GeomEntityConcept<FacetedVolume>, "FacetedVolume must satisfy GeomEntityConcept");

} // namespace gecko
