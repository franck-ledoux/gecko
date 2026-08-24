#pragma once
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <gecko/mesh/UnstructuredMesh.h>

namespace gecko::io {

    /**
     * @brief VTK legacy cell-type codes for the edges/faces/cells of a given @p TopologyTraits.
     * @tparam TopologyTraits gecko::SimplicialTraits or gecko::CubicTraits.
     */
    template<typename TopologyTraits>
    struct VtkElementTraits;

    /** @brief VTK legacy element types for a triangular/tetrahedral (SimplicialTraits) mesh. */
    template<>
    struct VtkElementTraits<SimplicialTraits> {
        /** @brief VTK_LINE. */
        static constexpr int edge_type = 3;
        /** @brief VTK_TRIANGLE. */
        static constexpr int face_type = 5;
        /** @brief VTK_TETRA. */
        static constexpr int cell_type = 10;
    };

    /** @brief VTK legacy element types for a quadrangular/hexahedral (CubicTraits) mesh. */
    template<>
    struct VtkElementTraits<CubicTraits> {
        /** @brief VTK_LINE. */
        static constexpr int edge_type = 3;
        /** @brief VTK_QUAD. */
        static constexpr int face_type = 9;
        /** @brief VTK_HEXAHEDRON. */
        static constexpr int cell_type = 12;
    };

    /**
     * @class VtkMeshWriter
     * @brief Writes a gecko::UnstructuredMesh to a VTK legacy ASCII file (`DATASET
     * UNSTRUCTURED_GRID`), mirroring gecko::io::GmshMeshWriter's shape for the VTK format.
     * @tparam TopologyTraits Topology of the mesh structure to write: gecko::SimplicialTraits
     * (faces written as VTK_TRIANGLE, cells as VTK_TETRA) or gecko::CubicTraits (faces written as
     * VTK_QUAD, cells as VTK_HEXAHEDRON) — see VtkElementTraits.
     */
    template<typename TopologyTraits>
    class VtkMeshWriter {
    public:
        /**
         * @brief Writes a gecko::UnstructuredMesh\<TopologyTraits\> to a VTK legacy ASCII file.
         *
         * Every node becomes a `POINTS` entry. Every edge/face/cell becomes one entry of the
         * dataset's single, unified `CELLS`/`CELL_TYPES` arrays — VTK legacy format has no notion of
         * grouping elements by topological dimension the way Gmsh's `$Elements` sections do, so
         * edges, faces and cells are all written into the same 2 arrays, in that order. A node
         * carrying no edge/face/cell of its own is still present in `POINTS` but contributes no
         * `CELLS` entry (VTK legacy has no standalone "point element").
         *
         * @p ANodeIntVariables optionally names `Variable<Int>` node variables already registered on
         * @p mesh (see `UnstructuredMesh::add_variable`) to write as `POINT_DATA`/`SCALARS` sections —
         * e.g. `Blocking::NODE_CLASSIFICATION_DIM_VARIABLE`/`NODE_CLASSIFICATION_TAG_VARIABLE`. A name
         * not actually present on @p mesh is silently skipped, so callers can pass a fixed list
         * without checking each one first. Nothing is written when the list is empty or none of its
         * names are present — this writer stays a general VTK legacy exporter, agnostic of what any
         * particular variable name means; it is the caller's job to know which variables exist.
         *
         * @p AElementIntVariables does the same for `CELL_DATA`, and has one wrinkle the point
         * side does not: VTK's cell data covers the single `CELLS` array, so one array has to span
         * edges, faces and cells together. A name is therefore looked up on all 3 registries and the
         * values concatenated in that same order, with -1 standing in wherever the name is not
         * registered for a dimension — `block_id` on an edge, say. A name present on none of the 3
         * is skipped, exactly as on the point side.
         *
         * @param path Path of the .vtk file to write (overwritten if it already exists).
         * @param mesh Mesh structure to write.
         * @param ANodeIntVariables Names of `Variable<Int>` node variables to write as `POINT_DATA`.
         * @param AElementIntVariables Names of `Variable<Int>` edge/face/cell variables to write as
         *        `CELL_DATA`.
         * @throw std::runtime_error if the file cannot be opened for writing.
         */
        static void write(const std::string &path,
                          const UnstructuredMesh<TopologyTraits> &mesh,
                          const std::vector<std::string> &ANodeIntVariables = {},
                          const std::vector<std::string> &AElementIntVariables = {}) {
            std::ofstream file(path);
            if (!file) {
                throw std::runtime_error("gecko::io::VtkMeshWriter::write: cannot open file '" + path +
                                         "' for writing");
            }
            file << std::setprecision(std::numeric_limits<Float>::max_digits10);

            file << "# vtk DataFile Version 3.0\n";
            file << "gecko blocking export\n";
            file << "ASCII\n";
            file << "DATASET UNSTRUCTURED_GRID\n";

            file << "POINTS " << mesh.nb_nodes() << " double\n";
            for (UInt i = 0; i < mesh.nb_nodes(); ++i) {
                const auto &p = mesh.node(NodeId{i});
                file << p.x() << " " << p.y() << " " << p.z() << "\n";
            }

            const UInt nb_elements = mesh.nb_edges() + mesh.nb_faces() + mesh.nb_cells();
            const UInt total_ints = mesh.nb_edges() * (1 + 2) + mesh.nb_faces() * (1 + TopologyTraits::Face::Nodes) +
                                    mesh.nb_cells() * (1 + TopologyTraits::Cell::Nodes);

            file << "CELLS " << nb_elements << " " << total_ints << "\n";
            for (UInt i = 0; i < mesh.nb_edges(); ++i) {
                file << 2;
                for (auto n : mesh.edge_nodes(EdgeId{i})) {
                    file << " " << n.value;
                }
                file << "\n";
            }
            for (UInt i = 0; i < mesh.nb_faces(); ++i) {
                file << TopologyTraits::Face::Nodes;
                for (auto n : mesh.face_nodes(FaceId{i})) {
                    file << " " << n.value;
                }
                file << "\n";
            }
            for (UInt i = 0; i < mesh.nb_cells(); ++i) {
                file << TopologyTraits::Cell::Nodes;
                for (auto n : mesh.cell_nodes(CellId{i})) {
                    file << " " << n.value;
                }
                file << "\n";
            }

            file << "CELL_TYPES " << nb_elements << "\n";
            for (UInt i = 0; i < mesh.nb_edges(); ++i) {
                file << VtkElementTraits<TopologyTraits>::edge_type << "\n";
            }
            for (UInt i = 0; i < mesh.nb_faces(); ++i) {
                file << VtkElementTraits<TopologyTraits>::face_type << "\n";
            }
            for (UInt i = 0; i < mesh.nb_cells(); ++i) {
                file << VtkElementTraits<TopologyTraits>::cell_type << "\n";
            }

            bool wrote_point_data_header = false;
            for (const std::string &name : ANodeIntVariables) {
                if (!mesh.template has_variable<CellType::Node>(name)) continue;
                if (!wrote_point_data_header) {
                    file << "POINT_DATA " << mesh.nb_nodes() << "\n";
                    wrote_point_data_header = true;
                }
                const auto &values = mesh.template get_variable<Int, CellType::Node>(name);
                file << "SCALARS " << name << " int 1\n";
                file << "LOOKUP_TABLE default\n";
                for (UInt i = 0; i < mesh.nb_nodes(); ++i) {
                    file << values[i] << "\n";
                }
            }

            bool wrote_cell_data_header = false;
            for (const std::string &name : AElementIntVariables) {
                const bool on_edges = mesh.template has_variable<CellType::Edge>(name);
                const bool on_faces = mesh.template has_variable<CellType::Face>(name);
                const bool on_cells = mesh.template has_variable<CellType::Cell>(name);
                if (!on_edges && !on_faces && !on_cells) continue;
                if (!wrote_cell_data_header) {
                    file << "CELL_DATA " << nb_elements << "\n";
                    wrote_cell_data_header = true;
                }
                file << "SCALARS " << name << " int 1\n";
                file << "LOOKUP_TABLE default\n";
                // In the order the CELLS array was written: edges, then faces, then cells.
                if (on_edges) {
                    const auto &values = mesh.template get_variable<Int, CellType::Edge>(name);
                    for (UInt i = 0; i < mesh.nb_edges(); ++i) {
                        file << values[i] << "\n";
                    }
                } else {
                    for (UInt i = 0; i < mesh.nb_edges(); ++i) {
                        file << "-1\n";
                    }
                }
                if (on_faces) {
                    const auto &values = mesh.template get_variable<Int, CellType::Face>(name);
                    for (UInt i = 0; i < mesh.nb_faces(); ++i) {
                        file << values[i] << "\n";
                    }
                } else {
                    for (UInt i = 0; i < mesh.nb_faces(); ++i) {
                        file << "-1\n";
                    }
                }
                if (on_cells) {
                    const auto &values = mesh.template get_variable<Int, CellType::Cell>(name);
                    for (UInt i = 0; i < mesh.nb_cells(); ++i) {
                        file << values[i] << "\n";
                    }
                } else {
                    for (UInt i = 0; i < mesh.nb_cells(); ++i) {
                        file << "-1\n";
                    }
                }
            }
        }
    };

    /** @brief Writes a triangular/tetrahedral mesh (see gecko::SimplicialMesh) to a VTK legacy ASCII file. */
    using SimplicialVtkWriter = VtkMeshWriter<SimplicialTraits>;
    /** @brief Writes a quadrangular/hexahedral mesh (see gecko::CubicMesh) to a VTK legacy ASCII file. */
    using CubicVtkWriter = VtkMeshWriter<CubicTraits>;

} // namespace gecko::io
