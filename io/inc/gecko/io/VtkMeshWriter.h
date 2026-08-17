#pragma once
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

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
         * `CELLS` entry (VTK legacy has no standalone "point element"). No group/tag data is written
         * (VTK legacy's optional `POINT_DATA`/`CELL_DATA` scalar sections are out of scope for this
         * export — see issue #22's "VTK legacy export of the block structure" ask, satisfied by
         * `write(path, blocking.to_mesh(1))`: 1 subdivision reproduces exactly the block corners,
         * one quad/hex per top-cell).
         *
         * @param path Path of the .vtk file to write (overwritten if it already exists).
         * @param mesh Mesh structure to write.
         * @throw std::runtime_error if the file cannot be opened for writing.
         */
        static void write(const std::string &path, const UnstructuredMesh<TopologyTraits> &mesh) {
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
        }
    };

    /** @brief Writes a triangular/tetrahedral mesh (see gecko::SimplicialMesh) to a VTK legacy ASCII file. */
    using SimplicialVtkWriter = VtkMeshWriter<SimplicialTraits>;
    /** @brief Writes a quadrangular/hexahedral mesh (see gecko::CubicMesh) to a VTK legacy ASCII file. */
    using CubicVtkWriter = VtkMeshWriter<CubicTraits>;

} // namespace gecko::io
