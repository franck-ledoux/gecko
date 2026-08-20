#pragma once

#include <algorithm>
#include <map>
#include <memory>
#include <ranges>
#include <span>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <gecko/geom/FacetedEntities.h>
#include <gecko/geom/FacetedEntityRefVariant.h>
#include <gecko/geom_itf/Concepts.h>
#include <gecko/io/GmshCommon.h>
#include <gecko/io/GmshMeshReader.h>
#include <gecko/mesh/UnstructuredMesh.h>
#include <gecko/utils/Groups.h>

namespace gecko {

    /**
     * @class FacetedGeometry
     * @brief Faceted (mesh-backed) implementation of GeomModelConcept: a CAD model reconstructed
     * from a Gmsh MSH file containing only triangles (a boundary representation: surfaces, curves
     * and vertices) or tetrahedra (which additionally reconstructs volumes).
     *
     * Every edge/face/cell sharing the same Gmsh elementary entity tag (see gecko::io::
     * ENTITY_TAG_VARIABLE) is grouped into one FacetedCurve/FacetedSurface/FacetedVolume;
     * FacetedVertex entities come from Gmsh point elements (type 15), read explicitly rather than
     * inferred topologically. Elements with entity tag 0 (no elementary entity) are ignored. The
     * model's groups() are exactly the file's `$PhysicalNames` (gecko::io::PHYSICAL_GROUP_VARIABLE)
     * — no synthetic per-entity-tag groups are registered, since GroupRegistry::add_group()
     * dedupes by (name, dimension) and a synthetic name could silently collide with an unrelated
     * physical group of the same name.
     */
    class FacetedGeometry {
    public:
        /**
         * @brief Constructor. Reads the backing mesh from a Gmsh MSH file and reconstructs every
         * geometric entity from its elements' elementary entity tags.
         *
         * Does **not** call `UnstructuredMesh::build_connectivity()`: that computes cell-to-cell
         * and face-to-face adjacency under a 2-manifold assumption (it throws if an edge borders
         * more than 2 faces), which a legitimate multi-volume B-Rep file routinely violates — an
         * edge where an internal partition surface between two volumes meets the exterior
         * boundary is shared by 3+ triangles from 3+ different surfaces, which is valid CAD
         * topology, not mesh corruption. FacetedGeometry never needs that adjacency anyway: every
         * entity is reconstructed purely from elementary entity tags (see build_entities()).
         *
         * @param path Path to the .msh file (triangles-only for a boundary rep, or tetrahedra to
         *        also get volumes).
         * @throw std::runtime_error if the file cannot be read (see GmshMeshReader::read()).
         */
        explicit FacetedGeometry(const std::string &path) {
            auto data = io::SimplicialMeshReader::read(path);
            m_mesh = std::make_unique<SimplicialMesh>(std::move(data.mesh));
            m_groups = std::move(data.groups);
            build_entities();
        }

        /** @brief Copy is disabled to avoid uncontrolled reallocations of the geometry. */
        FacetedGeometry(const FacetedGeometry &) = delete;
        /** @brief Copy is disabled to avoid uncontrolled reallocations of the geometry. @return Reference to this FacetedGeometry. */
        FacetedGeometry &operator=(const FacetedGeometry &) = delete;

        /** @brief Move constructor. */
        FacetedGeometry(FacetedGeometry &&) noexcept = default;
        /** @brief Move assignment. @return Reference to this FacetedGeometry. */
        FacetedGeometry &operator=(FacetedGeometry &&) noexcept = default;

        /** @brief Gets the number of vertices in the model. @return The vertex count. */
        [[nodiscard]] std::size_t nb_vertices() const noexcept { return m_vertices.size(); }
        /** @brief Gets the number of curves in the model. @return The curve count. */
        [[nodiscard]] std::size_t nb_curves() const noexcept { return m_curves.size(); }
        /** @brief Gets the number of surfaces in the model. @return The surface count. */
        [[nodiscard]] std::size_t nb_surfaces() const noexcept { return m_surfaces.size(); }
        /** @brief Gets the number of volumes in the model. @return The volume count. */
        [[nodiscard]] std::size_t nb_volumes() const noexcept { return m_volumes.size(); }

        /** @brief Gets a zero-copy view onto every vertex of the model. @return A span of vertices. */
        [[nodiscard]] std::span<const FacetedVertex> vertices() const noexcept { return m_vertices; }
        /** @brief Gets a zero-copy view onto every curve of the model. @return A span of curves. */
        [[nodiscard]] std::span<const FacetedCurve> curves() const noexcept { return m_curves; }
        /** @brief Gets a zero-copy view onto every surface of the model. @return A span of surfaces. */
        [[nodiscard]] std::span<const FacetedSurface> surfaces() const noexcept { return m_surfaces; }
        /** @brief Gets a zero-copy view onto every volume of the model. @return A span of volumes. */
        [[nodiscard]] std::span<const FacetedVolume> volumes() const noexcept { return m_volumes; }

        /**
         * @brief Finds the vertex reconstructed from a given Gmsh elementary entity tag.
         * @param tag Elementary entity tag to look up.
         * @return A pointer to the matching vertex, or nullptr if no vertex has that tag.
         */
        [[nodiscard]] const FacetedVertex *vertex_by_tag(Int tag) const {
            const auto it = m_vertex_by_tag.find(tag);
            return it != m_vertex_by_tag.end() ? it->second : nullptr;
        }

        /**
         * @brief Finds the curve reconstructed from a given Gmsh elementary entity tag.
         * @param tag Elementary entity tag to look up.
         * @return A pointer to the matching curve, or nullptr if no curve has that tag.
         */
        [[nodiscard]] const FacetedCurve *curve_by_tag(Int tag) const {
            const auto it = m_curve_by_tag.find(tag);
            return it != m_curve_by_tag.end() ? it->second : nullptr;
        }

        /**
         * @brief Finds the surface reconstructed from a given Gmsh elementary entity tag.
         * @param tag Elementary entity tag to look up.
         * @return A pointer to the matching surface, or nullptr if no surface has that tag.
         */
        [[nodiscard]] const FacetedSurface *surface_by_tag(Int tag) const {
            const auto it = m_surface_by_tag.find(tag);
            return it != m_surface_by_tag.end() ? it->second : nullptr;
        }

        /**
         * @brief Finds the volume reconstructed from a given Gmsh elementary entity tag.
         * @param tag Elementary entity tag to look up.
         * @return A pointer to the matching volume, or nullptr if no volume has that tag.
         */
        [[nodiscard]] const FacetedVolume *volume_by_tag(Int tag) const {
            const auto it = m_volume_by_tag.find(tag);
            return it != m_volume_by_tag.end() ? it->second : nullptr;
        }

        /**
         * @brief Gives access to the backing mesh every entity of this model was reconstructed
         * from — the triangles/tetrahedra themselves, as opposed to the entities grouping them.
         * Needed to display the model (see biy/) or otherwise consume its raw facets.
         * @return The backing mesh.
         */
        [[nodiscard]] const SimplicialMesh &mesh() const noexcept { return *m_mesh; }

        /**
         * @brief Gets a zero-copy view onto every physical group of the model, regardless of dimension.
         * @return A span of groups.
         */
        [[nodiscard]] std::span<const Group> groups() const noexcept { return m_groups.groups(); }

        /**
         * @brief Gets the physical groups of the model restricted to one dimension.
         * @param dim Dimension to filter groups on.
         * @return A lazy view over the matching groups.
         */
        [[nodiscard]] auto groups(GroupDim dim) const {
            return groups() | std::views::filter([dim](const Group &g) { return g.dimension == dim; });
        }

        /**
         * @brief Gets the (possibly dimension-mixed) entities belonging to a physical group.
         * @param gid Id of the group to query.
         * @return A span of type-erased entity references, empty if @p gid is invalid or unknown.
         */
        [[nodiscard]] std::span<const FacetedEntityRefVariant> entities(GroupId gid) const {
            if (!gid.is_valid() || gid.value >= m_entities_by_group.size()) {
                return {};
            }
            return m_entities_by_group[gid.value];
        }

        /**
         * @brief Gets the entity (@p dim, @p tag) itself plus every higher-dimensional entity of
         * the model containing it — a vertex's curves and surfaces, a curve's surfaces, and the
         * model's volumes.
         *
         * This is the model's B-Rep incidence, which Gmsh `$Elements` does not state anywhere:
         * every entity here is reconstructed independently from its elementary tag (see
         * build_entities()), so containment is instead recovered **exactly** from shared backing
         * mesh nodes — a vertex lies on a curve when its own mesh node is one of the curve's, and a
         * curve lies on a surface when every one of its mesh nodes is one of the surface's. Being
         * integer set membership, this involves no tolerance and cannot be ambiguous, unlike a
         * proximity test.
         *
         * Used by `Blocking::classify()` to infer what a block edge/face is classified on from the
         * classification of its own boundary — see its "lowest common containing entity" rule.
         *
         * @note Assumes the faceting is conformal: a curve's mesh nodes must be *the same node ids*
         *       as its bounding surface's. That holds for any `.msh` Gmsh produced from a single
         *       CAD model, but not for one assembled from independently meshed parts — where this
         *       simply reports fewer containments (never wrong ones), and callers fall back to
         *       their own geometric heuristics.
         *
         * @param dim Dimension of the entity to query.
         * @param tag `entity_tag()` of the entity to query.
         * @return The containing entities as (dimension, entity_tag) pairs, starting with
         *         (@p dim, @p tag) itself; empty if no entity of that dimension has that tag.
         */
        [[nodiscard]] std::span<const std::pair<GroupDim, Int>> containing_entities(GroupDim dim, Int tag) const {
            const auto it = m_containing_entities.find({dim, tag});
            return it != m_containing_entities.end() ? std::span<const std::pair<GroupDim, Int>>(it->second)
                                                     : std::span<const std::pair<GroupDim, Int>>{};
        }

    private:
        /**
         * @brief Registers a type-erased reference to @p entity into its physical group's entity
         * list, if it has a valid one.
         * @param gid Physical group of the entity (from its first constituent element), or an
         *        invalid GroupId if it has none.
         * @param entity Type-erased pointer to the entity to register.
         */
        void register_in_group(GroupId gid, FacetedEntityRefVariant entity) {
            if (gid.is_valid()) {
                m_entities_by_group[gid.value].push_back(entity);
            }
        }

        /**
         * @brief Groups the mesh's nodes/edges/faces/cells by their Gmsh elementary entity tag
         * (tag 0 excluded) into one FacetedVertex/Curve/Surface/Volume per distinct tag, then
         * indexes each by the physical group of its first constituent element.
         *
         * @note Assumes every element making up one geometric entity shares the same physical tag
         *       (standard Gmsh practice) and that a vertex's tag is carried by exactly one point
         *       element.
         */
        void build_entities() {
            std::map<Int, std::vector<NodeId>> nodes_by_tag;
            std::map<Int, std::vector<EdgeId>> edges_by_tag;
            std::map<Int, std::vector<FaceId>> faces_by_tag;
            std::map<Int, std::vector<CellId>> cells_by_tag;

            auto &node_entity = m_mesh->get_variable<Int, CellType::Node>(std::string(io::ENTITY_TAG_VARIABLE));
            auto &edge_entity = m_mesh->get_variable<Int, CellType::Edge>(std::string(io::ENTITY_TAG_VARIABLE));
            auto &face_entity = m_mesh->get_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
            auto &cell_entity = m_mesh->get_variable<Int, CellType::Cell>(std::string(io::ENTITY_TAG_VARIABLE));

            auto &node_group = m_mesh->get_variable<GroupId, CellType::Node>(std::string(io::PHYSICAL_GROUP_VARIABLE));
            auto &edge_group = m_mesh->get_variable<GroupId, CellType::Edge>(std::string(io::PHYSICAL_GROUP_VARIABLE));
            auto &face_group = m_mesh->get_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
            auto &cell_group = m_mesh->get_variable<GroupId, CellType::Cell>(std::string(io::PHYSICAL_GROUP_VARIABLE));

            for (UInt i = 0; i < m_mesh->nb_nodes(); ++i) {
                if (node_entity[i] != 0) {
                    nodes_by_tag[node_entity[i]].push_back(NodeId{i});
                }
            }
            for (UInt i = 0; i < m_mesh->nb_edges(); ++i) {
                if (edge_entity[i] != 0) {
                    edges_by_tag[edge_entity[i]].push_back(EdgeId{i});
                }
            }
            for (UInt i = 0; i < m_mesh->nb_faces(); ++i) {
                if (face_entity[i] != 0) {
                    faces_by_tag[face_entity[i]].push_back(FaceId{i});
                }
            }
            for (UInt i = 0; i < m_mesh->nb_cells(); ++i) {
                if (cell_entity[i] != 0) {
                    cells_by_tag[cell_entity[i]].push_back(CellId{i});
                }
            }

            m_entities_by_group.resize(m_groups.size());

            m_vertices.reserve(nodes_by_tag.size());
            for (auto &[tag, ids] : nodes_by_tag) {
                m_vertices.emplace_back(m_mesh.get(), ids.front(), tag);
                m_vertex_by_tag[tag] = &m_vertices.back();
                register_in_group(node_group[ids.front().value], &m_vertices.back());
            }

            m_curves.reserve(edges_by_tag.size());
            for (auto &[tag, ids] : edges_by_tag) {
                const auto gid = edge_group[ids.front().value];
                m_curves.emplace_back(m_mesh.get(), std::move(ids), tag);
                m_curve_by_tag[tag] = &m_curves.back();
                register_in_group(gid, &m_curves.back());
            }

            m_surfaces.reserve(faces_by_tag.size());
            for (auto &[tag, ids] : faces_by_tag) {
                const auto gid = face_group[ids.front().value];
                m_surfaces.emplace_back(m_mesh.get(), std::move(ids), tag);
                m_surface_by_tag[tag] = &m_surfaces.back();
                register_in_group(gid, &m_surfaces.back());
            }

            m_volumes.reserve(cells_by_tag.size());
            for (auto &[tag, ids] : cells_by_tag) {
                m_volumes.emplace_back(tag);
                m_volume_by_tag[tag] = &m_volumes.back();
                register_in_group(cell_group[ids.front().value], &m_volumes.back());
            }

            build_incidence();
        }

        /**
         * @brief Computes `containing_entities()`'s table from the entities' shared backing mesh
         * nodes — see that method for why incidence has to be recovered rather than read.
         */
        void build_incidence() {
            // Each curve's/surface's own node set, gathered once: containment is tested against
            // these repeatedly below.
            std::unordered_map<Int, std::unordered_set<UInt>> curve_nodes;
            std::unordered_map<Int, std::unordered_set<UInt>> surface_nodes;

            for (const auto &curve : m_curves) {
                auto &nodes = curve_nodes[curve.entity_tag()];
                for (const EdgeId e : curve.edges()) {
                    for (const NodeId n : m_mesh->edge_nodes(e)) {
                        nodes.insert(n.value);
                    }
                }
            }
            for (const auto &surface : m_surfaces) {
                auto &nodes = surface_nodes[surface.entity_tag()];
                for (const FaceId f : surface.faces()) {
                    for (const NodeId n : m_mesh->face_nodes(f)) {
                        nodes.insert(n.value);
                    }
                }
            }

            // Every entity contains itself, so a caller intersecting up-sets never gets an empty
            // result purely because two boundary cells share the same classification.
            const auto self_and = [this](GroupDim dim, Int tag) -> std::vector<std::pair<GroupDim, Int>> & {
                auto &up = m_containing_entities[{dim, tag}];
                up.emplace_back(dim, tag);
                return up;
            };

            for (const auto &vertex : m_vertices) {
                auto &up = self_and(GroupDim::Dim0, vertex.entity_tag());
                const UInt node = vertex.node().value;
                for (const auto &[tag, nodes] : curve_nodes) {
                    if (nodes.contains(node)) up.emplace_back(GroupDim::Dim1, tag);
                }
                for (const auto &[tag, nodes] : surface_nodes) {
                    if (nodes.contains(node)) up.emplace_back(GroupDim::Dim2, tag);
                }
                append_volumes(up);
            }

            for (const auto &curve : m_curves) {
                auto &up = self_and(GroupDim::Dim1, curve.entity_tag());
                const auto &own = curve_nodes[curve.entity_tag()];
                for (const auto &[tag, nodes] : surface_nodes) {
                    const bool inside = std::ranges::all_of(own, [&nodes](UInt n) { return nodes.contains(n); });
                    if (inside) up.emplace_back(GroupDim::Dim2, tag);
                }
                append_volumes(up);
            }

            for (const auto &surface : m_surfaces) {
                append_volumes(self_and(GroupDim::Dim2, surface.entity_tag()));
            }

            for (const auto &volume : m_volumes) {
                self_and(GroupDim::Dim3, volume.entity_tag());
            }
        }

        /**
         * @brief Appends every volume of the model to @p AUp. `FacetedVolume` is a stub covering
         * all of space (its `distance()` is always 0), so it contains every lower-dimensional
         * entity unconditionally — there is no node set to test against.
         * @param AUp The containing-entity list to append to.
         */
        void append_volumes(std::vector<std::pair<GroupDim, Int>> &AUp) const {
            for (const auto &volume : m_volumes) {
                AUp.emplace_back(GroupDim::Dim3, volume.entity_tag());
            }
        }

        std::unique_ptr<SimplicialMesh> m_mesh;
        GroupRegistry m_groups;
        std::vector<FacetedVertex> m_vertices;
        std::vector<FacetedCurve> m_curves;
        std::vector<FacetedSurface> m_surfaces;
        std::vector<FacetedVolume> m_volumes;
        std::unordered_map<Int, const FacetedVertex *> m_vertex_by_tag;
        std::unordered_map<Int, const FacetedCurve *> m_curve_by_tag;
        std::unordered_map<Int, const FacetedSurface *> m_surface_by_tag;
        std::unordered_map<Int, const FacetedVolume *> m_volume_by_tag;
        std::vector<std::vector<FacetedEntityRefVariant>> m_entities_by_group;
        /** @brief `containing_entities()`'s precomputed table. A `std::map` rather than a hash map
         * purely so `(GroupDim, Int)` works as a key without a hand-written hash — it is built once
         * and queried a handful of times per classified cell. */
        std::map<std::pair<GroupDim, Int>, std::vector<std::pair<GroupDim, Int>>> m_containing_entities;
    };
    static_assert(GeomModelConcept<FacetedGeometry>, "FacetedGeometry must satisfy GeomModelConcept");

} // namespace gecko
