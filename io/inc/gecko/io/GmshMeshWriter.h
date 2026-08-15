#pragma once
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <gecko/io/GmshCommon.h>
#include <gecko/mesh/UnstructuredMesh.h>
#include <gecko/utils/Groups.h>

namespace gecko::io {

    /**
     * @class GmshMeshWriter
     * @brief Writes a gecko::UnstructuredMesh (mesh structure) to a Gmsh MSH 2.2 ASCII file.
     * @tparam TopologyTraits Topology of the mesh structure to write: gecko::SimplicialTraits
     * (faces written as 3-node triangles, cells as 4-node tetrahedra) or gecko::CubicTraits
     * (faces written as 4-node quadrangles, cells as 8-node hexahedra).
     */
    template<typename TopologyTraits>
    class GmshMeshWriter {
    public:
        /**
         * @brief Writes a gecko::UnstructuredMesh\<TopologyTraits\> to a Gmsh MSH 2.2 ASCII file.
         *
         * Nodes and elements are written with fresh, contiguous 1-based ids (any original Gmsh
         * numbering is not preserved, since UnstructuredMesh does not keep it). Every node is
         * written to `$Nodes` regardless of tags, but only nodes carrying a physical group and/or
         * entity tag also get a 1-node point element (Gmsh type 15) in `$Elements` — untagged nodes
         * get no point element, since they are implied by the mesh's other elements. Edges are
         * written as 2-node lines (Gmsh type 1), faces and cells using
         * GmshElementTraits\<TopologyTraits\>::face_type / ::cell_type. For every point/edge/face/
         * cell, if @p mesh has a `Variable<GroupId>` named gecko::io::PHYSICAL_GROUP_VARIABLE its
         * group is written as the element's physical tag (element tag 0, using `GroupId::value + 1`
         * as the Gmsh physical tag number, so tags need not match those of a file this mesh was
         * originally read from), and if @p mesh has a `Variable<Int>` named
         * gecko::io::ENTITY_TAG_VARIABLE its value is written as the element's elementary entity tag
         * (element tag 1); when only one of the two is present, the same value is written for both
         * tags (see GmshCommon.h and GmshMeshReader::read()). @p groups is written to
         * `$PhysicalNames`.
         *
         * @param path Path of the .msh file to write (overwritten if it already exists).
         * @param mesh Mesh structure to write.
         * @param groups Registry of the mesh's physical groups (pass a default-constructed
         *        GroupRegistry{} if the mesh has none).
         * @throw std::runtime_error if the file cannot be opened for writing.
         */
        static void
        write(const std::string &path, const UnstructuredMesh<TopologyTraits> &mesh, const GroupRegistry &groups) {
            std::ofstream file(path);
            if (!file) {
                throw std::runtime_error("gecko::io::GmshMeshWriter::write: cannot open file '" + path +
                                         "' for writing");
            }
            file << std::setprecision(std::numeric_limits<Float>::max_digits10);

            file << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";

            if (groups.size() > 0) {
                file << "$PhysicalNames\n" << groups.size() << "\n";
                for (UInt i = 0; i < groups.size(); ++i) {
                    const auto &g = groups.get_group(GroupId{i});
                    file << static_cast<int>(g.dimension) << " " << (i + 1) << " \"" << g.name << "\"\n";
                }
                file << "$EndPhysicalNames\n";
            }

            file << "$Nodes\n" << mesh.nb_nodes() << "\n";
            for (UInt i = 0; i < mesh.nb_nodes(); ++i) {
                const auto &p = mesh.node(NodeId{i});
                file << (i + 1) << " " << p.x() << " " << p.y() << " " << p.z() << "\n";
            }
            file << "$EndNodes\n";

            const auto group_name = std::string(PHYSICAL_GROUP_VARIABLE);
            const auto entity_name = std::string(ENTITY_TAG_VARIABLE);

            const bool has_node_group = mesh.template has_variable<CellType::Node>(group_name);
            const bool has_node_entity = mesh.template has_variable<CellType::Node>(entity_name);
            std::vector<std::tuple<UInt, GroupId, Int>> tagged_nodes;
            if (has_node_group || has_node_entity) {
                for (UInt i = 0; i < mesh.nb_nodes(); ++i) {
                    const auto gid =
                        has_node_group ? mesh.template get_variable<GroupId, CellType::Node>(group_name)[i] : GroupId{};
                    const auto entity =
                        has_node_entity ? mesh.template get_variable<Int, CellType::Node>(entity_name)[i] : Int{0};
                    if (gid.is_valid() || entity != 0) {
                        tagged_nodes.emplace_back(i, gid, entity);
                    }
                }
            }

            file << "$Elements\n"
                 << (tagged_nodes.size() + mesh.nb_edges() + mesh.nb_faces() + mesh.nb_cells()) << "\n";
            UInt elm_number = 1;

            for (const auto &[i, gid, entity] : tagged_nodes) {
                file << elm_number << " " << POINT_TYPE;
                write_tags(file, gid, entity);
                file << " " << (i + 1) << "\n";
                ++elm_number;
            }

            const bool has_edge_group = mesh.template has_variable<CellType::Edge>(group_name);
            const bool has_edge_entity = mesh.template has_variable<CellType::Edge>(entity_name);
            for (UInt i = 0; i < mesh.nb_edges(); ++i, ++elm_number) {
                const auto gid =
                    has_edge_group ? mesh.template get_variable<GroupId, CellType::Edge>(group_name)[i] : GroupId{};
                const auto entity =
                    has_edge_entity ? mesh.template get_variable<Int, CellType::Edge>(entity_name)[i] : Int{0};
                file << elm_number << " " << EDGE_TYPE;
                write_tags(file, gid, entity);
                for (auto n : mesh.edge_nodes(EdgeId{i})) {
                    file << " " << (n.value + 1);
                }
                file << "\n";
            }

            const bool has_face_group = mesh.template has_variable<CellType::Face>(group_name);
            const bool has_face_entity = mesh.template has_variable<CellType::Face>(entity_name);
            for (UInt i = 0; i < mesh.nb_faces(); ++i, ++elm_number) {
                const auto gid =
                    has_face_group ? mesh.template get_variable<GroupId, CellType::Face>(group_name)[i] : GroupId{};
                const auto entity =
                    has_face_entity ? mesh.template get_variable<Int, CellType::Face>(entity_name)[i] : Int{0};
                file << elm_number << " " << GmshElementTraits<TopologyTraits>::face_type;
                write_tags(file, gid, entity);
                for (auto n : mesh.face_nodes(FaceId{i})) {
                    file << " " << (n.value + 1);
                }
                file << "\n";
            }

            const bool has_cell_group = mesh.template has_variable<CellType::Cell>(group_name);
            const bool has_cell_entity = mesh.template has_variable<CellType::Cell>(entity_name);
            for (UInt i = 0; i < mesh.nb_cells(); ++i, ++elm_number) {
                const auto gid =
                    has_cell_group ? mesh.template get_variable<GroupId, CellType::Cell>(group_name)[i] : GroupId{};
                const auto entity =
                    has_cell_entity ? mesh.template get_variable<Int, CellType::Cell>(entity_name)[i] : Int{0};
                file << elm_number << " " << GmshElementTraits<TopologyTraits>::cell_type;
                write_tags(file, gid, entity);
                for (auto n : mesh.cell_nodes(CellId{i})) {
                    file << " " << (n.value + 1);
                }
                file << "\n";
            }

            file << "$EndElements\n";
        }

    private:
        /**
         * @brief Writes an `$Elements` line's tag list: `0` if neither @p group nor @p entity is
         * set, otherwise `2 <physical> <entity>` (each falling back to the other when only one of
         * the two is available).
         * @param out Stream to write to.
         * @param group Physical group of the element, or an invalid GroupId if none.
         * @param entity Elementary entity tag of the element, or 0 if none.
         */
        static void write_tags(std::ostream &out, GroupId group, Int entity) {
            const bool has_physical = group.is_valid();
            const bool has_entity = entity != 0;
            if (!has_physical && !has_entity) {
                out << " 0";
                return;
            }
            const Int physical_tag = has_physical ? static_cast<Int>(group.value + 1) : entity;
            const Int entity_tag = has_entity ? entity : physical_tag;
            out << " 2 " << physical_tag << " " << entity_tag;
        }
    };

    /** @brief Writes a triangular/tetrahedral mesh (see gecko::SimplicialMesh) to a Gmsh MSH 2.2 ASCII file. */
    using SimplicialMeshWriter = GmshMeshWriter<SimplicialTraits>;
    /** @brief Writes a quadrangular/hexahedral mesh (see gecko::CubicMesh) to a Gmsh MSH 2.2 ASCII file. */
    using CubicMeshWriter = GmshMeshWriter<CubicTraits>;

} // namespace gecko::io
