#pragma once
#include <array>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>

#include <gecko/io/GmshCommon.h>
#include <gecko/mesh/UnstructuredMesh.h>
#include <gecko/utils/Groups.h>

namespace gecko::io {

    /**
     * @brief Result of reading a Gmsh .msh file: the mesh itself and the registry of its physical groups.
     * @tparam TopologyTraits Topology traits of the read mesh (SimplicialTraits or CubicTraits).
     */
    template<typename TopologyTraits>
    struct GmshMeshData {
        UnstructuredMesh<TopologyTraits> mesh; ///< The mesh structure read from the file.
        GroupRegistry groups;                  ///< The Gmsh physical groups (`$PhysicalNames`), indexed by GroupId.
    };

    /**
     * @class GmshMeshReader
     * @brief Reads a gecko::UnstructuredMesh (mesh structure) from a Gmsh MSH 2.2 ASCII file.
     * @tparam TopologyTraits Topology of the mesh structure to read: gecko::SimplicialTraits (reads
     * only 3-node triangles as faces and 4-node tetrahedra as cells) or gecko::CubicTraits (reads
     * only 4-node quadrangles as faces and 8-node hexahedra as cells).
     */
    template<typename TopologyTraits>
    class GmshMeshReader {
    public:
        /** @brief Topology (node count) of the faces this reader recognizes. */
        using Face = typename TopologyTraits::Face;
        /** @brief Topology (node count) of the cells this reader recognizes. */
        using Cell = typename TopologyTraits::Cell;

        /**
         * @brief Reads a gecko::UnstructuredMesh\<TopologyTraits\> from a Gmsh MSH 2.2 ASCII file.
         *
         * `$Elements` entries are read as follows: 1-node points (Gmsh type 15) set the tags of the
         * already-existing mesh Node they reference, 2-node lines (Gmsh type 1) become mesh Edges,
         * entries matching GmshElementTraits\<TopologyTraits\>::face_type become mesh Faces, and
         * entries matching GmshElementTraits\<TopologyTraits\>::cell_type become mesh Cells; every
         * other entry (elements of the other topology's types, e.g. quads found while reading with
         * TopologyTraits = SimplicialTraits) is ignored. For every point/edge/face/cell read, its
         * first element tag (the physical tag), when present, is interpreted as its Gmsh physical
         * group and stored in a `Variable<GroupId>` named gecko::io::PHYSICAL_GROUP_VARIABLE, and
         * its second element tag (the elementary entity tag), when present, is stored verbatim in a
         * `Variable<Int>` named gecko::io::ENTITY_TAG_VARIABLE (see GmshCommon.h). The physical
         * groups themselves (`$PhysicalNames`) are returned via GmshMeshData::groups. Cell/face/edge
         * adjacency (`UnstructuredMesh::build_connectivity()`) is **not** computed automatically;
         * call it explicitly on the returned mesh if needed.
         *
         * @param path Path to the .msh file to read.
         * @return The mesh structure and its physical group registry.
         * @throw std::runtime_error if the file cannot be opened, is not a supported MSH 2.x file,
         *        or a required section is malformed.
         */
        static GmshMeshData<TopologyTraits> read(const std::string &path) {
            std::ifstream file(path);
            if (!file) {
                throw std::runtime_error("gecko::io::GmshMeshReader::read: cannot open file '" + path + "'");
            }

            GmshMeshData<TopologyTraits> data;
            auto &node_group =
                data.mesh.template add_variable<GroupId, CellType::Node>(std::string(PHYSICAL_GROUP_VARIABLE));
            auto &edge_group =
                data.mesh.template add_variable<GroupId, CellType::Edge>(std::string(PHYSICAL_GROUP_VARIABLE));
            auto &face_group =
                data.mesh.template add_variable<GroupId, CellType::Face>(std::string(PHYSICAL_GROUP_VARIABLE));
            auto &cell_group =
                data.mesh.template add_variable<GroupId, CellType::Cell>(std::string(PHYSICAL_GROUP_VARIABLE));
            auto &node_entity = data.mesh.template add_variable<Int, CellType::Node>(std::string(ENTITY_TAG_VARIABLE));
            auto &edge_entity = data.mesh.template add_variable<Int, CellType::Edge>(std::string(ENTITY_TAG_VARIABLE));
            auto &face_entity = data.mesh.template add_variable<Int, CellType::Face>(std::string(ENTITY_TAG_VARIABLE));
            auto &cell_entity = data.mesh.template add_variable<Int, CellType::Cell>(std::string(ENTITY_TAG_VARIABLE));

            std::unordered_map<Int, NodeId> node_id_of;       // Gmsh node-number -> NodeId
            std::unordered_map<Int, GroupId> group_id_of_tag; // Gmsh physical tag -> GroupId

            const auto next_line = [&file](const char *section) -> std::string {
                std::string line;
                if (!std::getline(file, line)) {
                    throw std::runtime_error(std::string("gecko::io::GmshMeshReader::read: truncated ") + section +
                                             " section");
                }
                strip_cr(line);
                return line;
            };

            bool has_mesh_format = false;
            std::string line;
            while (std::getline(file, line)) {
                strip_cr(line);
                if (line == "$MeshFormat") {
                    std::istringstream iss(next_line("$MeshFormat"));
                    std::string version;
                    iss >> version;
                    if (version.empty() || version[0] != '2') {
                        throw std::runtime_error("gecko::io::GmshMeshReader::read: unsupported MSH version '" +
                                                 version + "', only MSH 2.x is supported");
                    }
                    has_mesh_format = true;
                } else if (line == "$PhysicalNames") {
                    const auto count = std::stoul(next_line("$PhysicalNames"));
                    for (std::size_t i = 0; i < count; ++i) {
                        std::istringstream iss(next_line("$PhysicalNames"));
                        int dim = 0;
                        Int tag = 0;
                        std::string name;
                        iss >> dim >> tag >> name;
                        const auto group_id = data.groups.add_group(unquote(name), to_group_dim(dim));
                        group_id_of_tag[tag] = group_id;
                    }
                } else if (line == "$Nodes") {
                    const auto count = std::stoul(next_line("$Nodes"));
                    for (std::size_t i = 0; i < count; ++i) {
                        std::istringstream iss(next_line("$Nodes"));
                        Int number = 0;
                        Float x = 0.0, y = 0.0, z = 0.0;
                        iss >> number >> x >> y >> z;
                        node_id_of[number] = data.mesh.add_node(x, y, z);
                    }
                } else if (line == "$Elements") {
                    const auto count = std::stoul(next_line("$Elements"));
                    for (std::size_t i = 0; i < count; ++i) {
                        std::istringstream iss(next_line("$Elements"));
                        Int elm_number = 0, elm_type = 0, num_tags = 0;
                        iss >> elm_number >> elm_type >> num_tags;

                        // Element tag 0 is the physical tag, tag 1 is the elementary entity tag;
                        // any further tags (e.g. partition ids) are read but not kept.
                        Int physical_tag = 0;
                        Int entity_tag = 0;
                        for (Int t = 0; t < num_tags; ++t) {
                            Int tag = 0;
                            iss >> tag;
                            if (t == 0) {
                                physical_tag = tag;
                            } else if (t == 1) {
                                entity_tag = tag;
                            }
                        }

                        GroupId group{};
                        if (num_tags > 0) {
                            if (auto it = group_id_of_tag.find(physical_tag); it != group_id_of_tag.end()) {
                                group = it->second;
                            }
                        }

                        if (elm_type == POINT_TYPE) {
                            const auto ids = read_node_ids<1>(iss, node_id_of);
                            node_group[ids[0].value] = group;
                            node_entity[ids[0].value] = entity_tag;
                        } else if (elm_type == EDGE_TYPE) {
                            const auto ids = read_node_ids<2>(iss, node_id_of);
                            auto e = data.mesh.add_edge(ids[0], ids[1]);
                            edge_group[e.value] = group;
                            edge_entity[e.value] = entity_tag;
                        } else if (elm_type == GmshElementTraits<TopologyTraits>::face_type) {
                            const auto ids = read_node_ids<Face::Nodes>(iss, node_id_of);
                            auto f = std::apply([&](auto... n) { return data.mesh.add_face(n...); }, ids);
                            face_group[f.value] = group;
                            face_entity[f.value] = entity_tag;
                        } else if (elm_type == GmshElementTraits<TopologyTraits>::cell_type) {
                            const auto ids = read_node_ids<Cell::Nodes>(iss, node_id_of);
                            auto c = std::apply([&](auto... n) { return data.mesh.add_cell(n...); }, ids);
                            cell_group[c.value] = group;
                            cell_entity[c.value] = entity_tag;
                        }
                        // Every other element type (points, or the other topology's face/cell types) is ignored.
                    }
                }
            }

            if (!has_mesh_format) {
                throw std::runtime_error("gecko::io::GmshMeshReader::read: missing $MeshFormat section");
            }

            return data;
        }

    private:
        /**
         * @brief Reads @p N whitespace-separated Gmsh node numbers from @p iss and resolves them to NodeId.
         * @tparam N Number of node ids to read.
         * @param iss Stream positioned right before the node numbers on an `$Elements` line.
         * @param node_id_of Map from Gmsh node-number to NodeId, built while reading `$Nodes`.
         * @return The @p N resolved node ids, in the order they appear in the file.
         */
        template<std::size_t N>
        static std::array<NodeId, N> read_node_ids(std::istringstream &iss,
                                                   const std::unordered_map<Int, NodeId> &node_id_of) {
            std::array<NodeId, N> ids{};
            for (std::size_t k = 0; k < N; ++k) {
                Int number = 0;
                iss >> number;
                ids[k] = node_id_of.at(number);
            }
            return ids;
        }

        /**
         * @brief Strips a trailing carriage-return character left in a line by a file with Windows line endings.
         * @param line Line to strip in place.
         */
        static void strip_cr(std::string &line) {
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
        }

        /**
         * @brief Converts a Gmsh `$PhysicalNames` dimension (0 to 3) to a GroupDim.
         * @param dim Gmsh dimension value.
         * @return The corresponding GroupDim, or GroupDim::Undefined if @p dim is not in [0, 3].
         */
        static GroupDim to_group_dim(int dim) {
            switch (dim) {
                case 0:
                    return GroupDim::Dim0;
                case 1:
                    return GroupDim::Dim1;
                case 2:
                    return GroupDim::Dim2;
                case 3:
                    return GroupDim::Dim3;
                default:
                    return GroupDim::Undefined;
            }
        }

        /**
         * @brief Strips the surrounding double quotes from a Gmsh-formatted physical name.
         * @param s String to unquote (e.g. `"Inlet"`).
         * @return @p s without its surrounding quotes (e.g. `Inlet`), or @p s unchanged if it isn't quoted.
         */
        static std::string unquote(const std::string &s) {
            if (s.size() >= 2 && s.front() == '"' && s.back() == '"') {
                return s.substr(1, s.size() - 2);
            }
            return s;
        }
    };

    /** @brief Reads a triangular/tetrahedral mesh (see gecko::SimplicialMesh) from a Gmsh MSH 2.2 ASCII file. */
    using SimplicialMeshReader = GmshMeshReader<SimplicialTraits>;
    /** @brief Reads a quadrangular/hexahedral mesh (see gecko::CubicMesh) from a Gmsh MSH 2.2 ASCII file. */
    using CubicMeshReader = GmshMeshReader<CubicTraits>;

} // namespace gecko::io
