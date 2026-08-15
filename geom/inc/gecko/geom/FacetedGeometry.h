#pragma once

#include <map>
#include <memory>
#include <ranges>
#include <span>
#include <string>
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
         * @param path Path to the .msh file (triangles-only for a boundary rep, or tetrahedra to
         *        also get volumes).
         * @throw std::runtime_error if the file cannot be read (see GmshMeshReader::read()).
         */
        explicit FacetedGeometry(const std::string &path) {
            auto data = io::SimplicialMeshReader::read(path);
            m_mesh = std::make_unique<SimplicialMesh>(std::move(data.mesh));
            m_groups = std::move(data.groups);
            m_mesh->build_connectivity();
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
                register_in_group(node_group[ids.front().value], &m_vertices.back());
            }

            m_curves.reserve(edges_by_tag.size());
            for (auto &[tag, ids] : edges_by_tag) {
                const auto gid = edge_group[ids.front().value];
                m_curves.emplace_back(m_mesh.get(), std::move(ids), tag);
                register_in_group(gid, &m_curves.back());
            }

            m_surfaces.reserve(faces_by_tag.size());
            for (auto &[tag, ids] : faces_by_tag) {
                const auto gid = face_group[ids.front().value];
                m_surfaces.emplace_back(m_mesh.get(), std::move(ids), tag);
                register_in_group(gid, &m_surfaces.back());
            }

            m_volumes.reserve(cells_by_tag.size());
            for (auto &[tag, ids] : cells_by_tag) {
                m_volumes.emplace_back(tag);
                register_in_group(cell_group[ids.front().value], &m_volumes.back());
            }
        }

        std::unique_ptr<SimplicialMesh> m_mesh;
        GroupRegistry m_groups;
        std::vector<FacetedVertex> m_vertices;
        std::vector<FacetedCurve> m_curves;
        std::vector<FacetedSurface> m_surfaces;
        std::vector<FacetedVolume> m_volumes;
        std::vector<std::vector<FacetedEntityRefVariant>> m_entities_by_group;
    };
    static_assert(GeomModelConcept<FacetedGeometry>, "FacetedGeometry must satisfy GeomModelConcept");

} // namespace gecko
