#pragma once

#include <cassert>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Triangle_3.h>
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
     * rather than copied coordinates. Closest-point/distance queries are accelerated by a
     * CGAL::AABB_tree built once at construction time from those triangles (same
     * Simple_cartesian<double> kernel + AABB_triangle_primitive pattern as gmds_core's
     * FACSurface), rather than a brute-force scan; CGAL's AABB_tree package is
     * GPL-3.0-or-later, combinable with Gecko's AGPL-3.0 license under GPLv3/AGPLv3 section 13.
     *
     * Move-only (matching CGAL::AABB_tree, which disables copy): the tree's primitives hold
     * iterators into #m_cgal_triangles, so this class must never be copied, and relies on
     * std::vector's move constructor leaving the underlying buffer address (and hence those
     * iterators) unchanged when moved.
     */
    class FacetedSurface {
    public:
        /**
         * @brief Constructor. Builds the AABB tree from @p faces' triangles.
         * @param mesh Backing mesh; must outlive this surface.
         * @param faces Ids of the mesh faces this surface is built from (must be non-empty).
         * @param entity_tag Gmsh elementary entity tag this surface was reconstructed from.
         */
        FacetedSurface(const SimplicialMesh *mesh, std::vector<FaceId> faces, Int entity_tag)
            : m_mesh(mesh),
              m_faces(std::move(faces)),
              m_entity_tag(entity_tag),
              m_cgal_triangles(build_cgal_triangles(m_mesh, m_faces)),
              m_aabb_tree(m_cgal_triangles.begin(), m_cgal_triangles.end()) {
            assert(!m_faces.empty() && "FacetedSurface: must be built from at least one face");
            m_aabb_tree.accelerate_distance_queries();
        }

        /**
         * @brief Gets the closest point of this surface to a query point.
         * @param p Query point.
         * @return The closest point of the surface (over every constituent triangle) to @p p.
         */
        [[nodiscard]] Point3d closest_point(const Point3d &p) const {
            const auto cp = m_aabb_tree.closest_point(CgalPoint(p.x(), p.y(), p.z()));
            return Point3d(cp.x(), cp.y(), cp.z());
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
        using CgalKernel = CGAL::Simple_cartesian<double>;
        using CgalPoint = CgalKernel::Point_3;
        using CgalTriangle = CgalKernel::Triangle_3;
        using CgalTrianglePrimitive = CGAL::AABB_triangle_primitive<CgalKernel, std::vector<CgalTriangle>::iterator>;
        using CgalTraits = CGAL::AABB_traits<CgalKernel, CgalTrianglePrimitive>;
        using CgalTree = CGAL::AABB_tree<CgalTraits>;

        /**
         * @brief Builds the list of CGAL triangles making up @p faces, for the AABB tree.
         * @param mesh Backing mesh to read node positions from.
         * @param faces Ids of the mesh faces to convert.
         * @return One CgalTriangle per face, in the same order.
         */
        static std::vector<CgalTriangle> build_cgal_triangles(const SimplicialMesh *mesh,
                                                              std::span<const FaceId> faces) {
            std::vector<CgalTriangle> triangles;
            triangles.reserve(faces.size());
            for (FaceId f : faces) {
                const auto nodes = mesh->face_nodes(f);
                const auto &a = mesh->node(nodes[0]);
                const auto &b = mesh->node(nodes[1]);
                const auto &c = mesh->node(nodes[2]);
                triangles.emplace_back(
                    CgalPoint(a.x(), a.y(), a.z()), CgalPoint(b.x(), b.y(), b.z()), CgalPoint(c.x(), c.y(), c.z()));
            }
            return triangles;
        }

        const SimplicialMesh *m_mesh;
        std::vector<FaceId> m_faces;
        Int m_entity_tag;
        std::vector<CgalTriangle> m_cgal_triangles;
        CgalTree m_aabb_tree;
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
