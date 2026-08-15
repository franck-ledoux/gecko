/*----------------------------------------------------------------------------*/
#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <limits>
#include <map>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
// -----------------------------------------------------------------------------
#include <gecko/math/Point3d.h>
#include <gecko/mesh/VariableRegistry.h>
#include <gecko/utils/Types.h>
// -----------------------------------------------------------------------------
// Traits to define mesh topology
// -----------------------------------------------------------------------------
namespace gecko {
    /** @brief Topology traits describing simplicial cells (tetrahedra / triangles). */
    struct SimplicialTraits {
        /** @brief Topology of a 3D cell (tetrahedron): 4 nodes, 4 faces, 6 edges. */
        struct Cell {
            static constexpr SizeT Nodes = 4; ///< Number of nodes per cell.
            static constexpr SizeT Faces = 4; ///< Number of faces per cell.
            static constexpr SizeT Edges = 6; ///< Number of edges per cell.
            /**
             * @brief Local node indices forming each of the 4 triangular faces of a tetrahedron.
             * Face i is the complement of local vertex i (a triangle uses the 3 other vertices);
             * this decomposition is combinatorially forced for any 4-node tetrahedron.
             */
            static constexpr std::array<std::array<SizeT, 3>, 4> FaceNodes = {
                {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}}};
            /**
             * @brief Local node indices forming each of the 6 edges of a tetrahedron: every pair of
             * its 4 vertices. This decomposition is combinatorially forced for any 4-node
             * tetrahedron (C(4,2) = 6).
             */
            static constexpr std::array<std::array<SizeT, 2>, 6> EdgeNodes = {
                {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
        };

        /** @brief Topology of a 2D face (triangle): 3 nodes, 3 edges. */
        struct Face {
            static constexpr SizeT Nodes = 3; ///< Number of nodes per face.
            static constexpr SizeT Edges = 3; ///< Number of edges per face.
            /** @brief Local node indices forming each of the 3 edges (sequential perimeter). */
            static constexpr std::array<std::array<SizeT, 2>, 3> EdgeNodes = {{{0, 1}, {1, 2}, {2, 0}}};
        };
    };
    static_assert(SimplicialTraits::Cell::FaceNodes[0].size() == SimplicialTraits::Face::Nodes,
                  "SimplicialTraits::Cell::FaceNodes rows must have Face::Nodes entries");
    static_assert(SimplicialTraits::Cell::EdgeNodes[0].size() == 2,
                  "SimplicialTraits::Cell::EdgeNodes rows must have 2 entries (an edge always has 2 nodes)");
    static_assert(SimplicialTraits::Face::EdgeNodes[0].size() == 2,
                  "SimplicialTraits::Face::EdgeNodes rows must have 2 entries (an edge always has 2 nodes)");

    // -----------------------------------------------------------------------------
    /** @brief Topology traits describing cubic cells (hexahedra / quadrilaterals). */
    struct CubicTraits {
        /** @brief Topology of a 3D cell (hexahedron): 8 nodes, 6 faces, 12 edges. */
        struct Cell {
            static constexpr SizeT Nodes = 8;  ///< Number of nodes per cell.
            static constexpr SizeT Faces = 6;  ///< Number of faces per cell.
            static constexpr SizeT Edges = 12; ///< Number of edges per cell.
            /**
             * @brief Local node indices forming each of the 6 quadrilateral faces of a hexahedron.
             * Follows the standard Gmsh/VTK HEX8 node ordering (nodes 0-3 form the bottom face in
             * perimeter order, nodes 4-7 form the top face directly above 0-3), with each face
             * wound for an outward-facing normal.
             */
            static constexpr std::array<std::array<SizeT, 4>, 6> FaceNodes = {{
                {0, 3, 2, 1}, // bottom
                {4, 5, 6, 7}, // top
                {0, 1, 5, 4}, // front
                {1, 2, 6, 5}, // right
                {2, 3, 7, 6}, // back
                {3, 0, 4, 7}  // left
            }};
            /**
             * @brief Local node indices forming each of the 12 edges of a hexahedron, using the
             * same HEX8 node ordering as #FaceNodes: 4 bottom perimeter edges, 4 top perimeter
             * edges, and 4 vertical edges connecting the bottom face to the top face.
             */
            static constexpr std::array<std::array<SizeT, 2>, 12> EdgeNodes = {{
                {0, 1},
                {1, 2},
                {2, 3},
                {3, 0}, // bottom perimeter
                {4, 5},
                {5, 6},
                {6, 7},
                {7, 4}, // top perimeter
                {0, 4},
                {1, 5},
                {2, 6},
                {3, 7} // vertical
            }};
        };

        /** @brief Topology of a 2D face (quadrilateral): 4 nodes, 4 edges. */
        struct Face {
            static constexpr SizeT Nodes = 4; ///< Number of nodes per face.
            static constexpr SizeT Edges = 4; ///< Number of edges per face.
            /** @brief Local node indices forming each of the 4 edges (sequential perimeter). */
            static constexpr std::array<std::array<SizeT, 2>, 4> EdgeNodes = {{{0, 1}, {1, 2}, {2, 3}, {3, 0}}};
        };
    };
    static_assert(CubicTraits::Cell::FaceNodes[0].size() == CubicTraits::Face::Nodes,
                  "CubicTraits::Cell::FaceNodes rows must have Face::Nodes entries");
    static_assert(CubicTraits::Cell::EdgeNodes[0].size() == 2,
                  "CubicTraits::Cell::EdgeNodes rows must have 2 entries (an edge always has 2 nodes)");
    static_assert(CubicTraits::Face::EdgeNodes[0].size() == 2,
                  "CubicTraits::Face::EdgeNodes rows must have 2 entries (an edge always has 2 nodes)");
    /*----------------------------------------------------------------------------*/
    /** @brief Dimension of a mesh cell type. */
    enum class CellType { Node, Edge, Face, Cell };

    /*----------------------------------------------------------------------------*/
    /**
     * @class NeighborRef
     * @brief Encoded reference to a mesh entity's neighbor across one of its local facets, plus
     * the signed-code convention used to store it compactly in m_c2c_raw/m_f2f_raw (see
     * UnstructuredMesh::build_connectivity()).
     *
     * Reused at two dimension pairs: for a Cell's neighbor across a local face (Cell3D = another
     * cell, Boundary2D = a Face materialized at that boundary facet), and, one dimension down, for
     * a Face's neighbor across a local edge (Cell3D = another face, Boundary2D = an Edge
     * materialized at that free/boundary edge). None means no same-dimension neighbor exists and
     * no lower-dimension entity was materialized for that boundary facet either.
     */
    struct NeighborRef {
        /** @brief Kind of entity a NeighborRef points to. */
        enum class Type { Cell3D, Boundary2D, None };

        Type type; ///< Whether the referenced entity is same-dimension, lower-dimension, or absent.
        Int index; ///< Index of the referenced entity within its own container; meaningless when type == None.

        /**
         * @brief Sentinel signed neighbor code meaning "no same-dimension neighbor, and no
         * lower-dimension entity materialized for this boundary facet either". Used both as the
         * placeholder value pushed by add_face()/add_cell() before build_connectivity() has run,
         * and as the final value left by build_connectivity() for such a facet.
         */
        static constexpr Int NO_NEIGHBOR = std::numeric_limits<Int>::min();

        /**
         * @brief Encodes a boundary entity index into the signed neighbor code convention.
         * @param boundary_idx Index of the boundary entity (>= 0).
         * @return A negative code identifying the boundary entity, decodable by decode().
         * @pre @p boundary_idx must be non-negative and strictly less than the maximum Int value
         *      (the encoding collides with #NO_NEIGHBOR at the very last representable value),
         *      checked via assert in Debug builds.
         */
        static Int encode_boundary(Int boundary_idx) {
            assert(boundary_idx >= 0 && boundary_idx < std::numeric_limits<Int>::max() &&
                   "NeighborRef::encode_boundary: index out of representable range");
            return -boundary_idx - 1;
        }

        /**
         * @brief Decodes a signed neighbor code into a NeighborRef.
         *
         * #NO_NEIGHBOR decodes to Type::None; otherwise a non-negative code refers to a
         * same-dimension neighbor index, and a negative code refers to an encoded lower-dimension
         * entity index (see encode_boundary()).
         *
         * @param code Signed neighbor code to decode.
         * @return The decoded NeighborRef, tagged as Cell3D, Boundary2D or None.
         */
        static NeighborRef decode(Int code) {
            if (code == NO_NEIGHBOR) {
                return {Type::None, 0};
            } else if (code >= 0) {
                return {Type::Cell3D, code};
            } else {
                return {Type::Boundary2D, -code - 1};
            }
        }
    };

    /**
     * @class UnstructuredMesh
     * @brief Unstructured mesh storing nodes, edges, faces and 3D cells whose topology
     * (number of nodes/faces per cell, nodes/edges per face) is given by @p TopologyTraits.
     * @tparam TopologyTraits Topology traits type (e.g. SimplicialTraits or CubicTraits)
     * describing the Cell and Face topologies.
     */
    /*----------------------------------------------------------------------------*/
    template<typename TopologyTraits>
    class UnstructuredMesh {
    public:
        /** @brief Topology (nodes/faces count) of a 3D cell, as given by @p TopologyTraits. */
        using Cell = typename TopologyTraits::Cell;
        /** @brief Topology (nodes/edges count) of a 2D face, as given by @p TopologyTraits. */
        using Face = typename TopologyTraits::Face;
        /**
         * @brief Constructor
         * @param capacity Initial capacity of the node container.
         */
        explicit UnstructuredMesh(SizeT capacity = 64) {
            //nodes data
            m_nodes.reserve(capacity);
            //edges data
            m_e2n.reserve(capacity);
            //faces data
            m_f2n.reserve(capacity);
            m_f2c_raw.reserve(capacity);
            m_f2f_raw.reserve(capacity);
            //cells data
            m_c2n.reserve(capacity);
            m_c2c_raw.reserve(capacity);
        }

        /** @brief Gets the number of nodes stored in the mesh. @return The node count. */
        [[nodiscard]] UInt nb_nodes() const { return m_nodes.size(); }
        /** @brief Gets the number of edges stored in the mesh. @return The edge count. */
        [[nodiscard]] UInt nb_edges() const { return m_e2n.size() / 2; }
        /** @brief Gets the number of faces stored in the mesh. @return The face count. */
        [[nodiscard]] UInt nb_faces() const { return m_f2n.size() / Face::Nodes; }
        /** @brief Gets the number of 3D cells stored in the mesh. @return The cell count. */
        [[nodiscard]] UInt nb_cells() const { return m_c2n.size() / Cell::Nodes; }

        /**
         * @brief Adds a new node to the mesh at the given position.
         * @param pnt Position of the new node.
         * @return The id of the newly created node.
         */
        NodeId add_node(const Point3d &pnt) {
            const auto index = m_nodes.size();
            m_nodes.push_back(pnt);
            m_node_variables.resize_all(m_nodes.size());
            return NodeId{static_cast<NodeId>(index)};
        }

        /**
         * @brief Adds a new node to the mesh from its coordinates.
         * @param x X coordinate of the new node.
         * @param y Y coordinate of the new node.
         * @param z Z coordinate of the new node (default: 0.0).
         * @return The id of the newly created node.
         */
        NodeId add_node(Float x, Float y, Float z = 0.0) { return add_node(Point3d{x, y, z}); }

        /**
         * @brief Adds a new edge to the mesh from its two node ids.
         * @param n0 First node id of the edge.
         * @param n1 Second node id of the edge.
         * @return The id of the newly created edge.
         * @pre @p n0 and @p n1 must be valid and refer to nodes already present in the mesh
         *      (id.value < nb_nodes()), checked via assert in Debug builds.
         */
        EdgeId add_edge(NodeId n0, NodeId n1) {
            assert(n0.is_valid() && n0.value < m_nodes.size() && "UnstructuredMesh::add_edge: node id out of range");
            assert(n1.is_valid() && n1.value < m_nodes.size() && "UnstructuredMesh::add_edge: node id out of range");

            const auto index = nb_edges();
            m_e2n.push_back(n0);
            m_e2n.push_back(n1);

            m_edge_variables.resize_all(nb_edges());
            return EdgeId{static_cast<EdgeId>(index)};
        }

        /**
         * @brief Adds a new face to the mesh from its node ids.
         * @tparam Args Pack of NodeId, exactly Face::Nodes of them.
         * @param nodes The Face::Nodes node ids defining the face, in order.
         * @return The id of the newly created face.
         * @pre Every node id in @p nodes must be valid and refer to a node already
         *      present in the mesh (id.value < nb_nodes()), checked via assert in Debug builds.
         * @note Enabled only if the number of arguments matches Face::Nodes.
         */
        template<typename... Args>
        FaceId add_face(Args... nodes) {
            // 1. Compile-time checks: right arity, and every argument is a NodeId
            static_assert(sizeof...(Args) == Face::Nodes,
                          "add_face(): the number of nodes does not match the face topology!");
            static_assert((std::is_same_v<Args, NodeId> && ...), "add_face(): all the arguments must be NodeId!");

            // 2. Append the Face::Nodes node ids to the flat node-to-face connectivity array
            const auto index = nb_faces();
            for (NodeId n : {nodes...}) {
                assert(n.is_valid() && n.value < m_nodes.size() && "UnstructuredMesh::add_face: node id out of range");
                m_f2n.push_back(n);
            }

            // 3. No neighbor is known yet across the new face's edges; real adjacency (matching
            // shared edges between faces) will be computed separately, see build_connectivity()
            for (SizeT i = 0; i < Face::Edges; ++i) {
                m_f2f_raw.push_back(NeighborRef::NO_NEIGHBOR);
            }
            // idem for faces
            m_f2c_raw.push_back(CellId::invalid_id);
            m_f2c_raw.push_back(CellId::invalid_id);

            m_face_variables.resize_all(nb_faces());
            return FaceId{static_cast<FaceId>(index)};
        }

        /**
         * @brief Adds a new 3D cell to the mesh from its node ids.
         * @tparam Args Pack of NodeId, exactly Cell::Nodes of them.
         * @param nodes The Cell::Nodes node ids defining the cell, in order.
         * @return The id of the newly created cell.
         * @pre Every node id in @p nodes must be valid and refer to a node already
         *      present in the mesh (id.value < nb_nodes()), checked via assert in Debug builds.
         * @note Enabled only if the number of arguments matches Cell::Nodes.
         */
        template<typename... Args>
        CellId add_cell(Args... nodes) {
            // 1. Compile-time checks: right arity, and every argument is a NodeId
            static_assert(sizeof...(Args) == Cell::Nodes,
                          "add_cell(): the number of nodes does not match the cell topology!");
            static_assert((std::is_same_v<Args, NodeId> && ...), "add_cell(): all the arguments must be NodeId!");

            // 2. Append the Cell::Nodes node ids to the flat node-to-cell connectivity array
            const auto index = nb_cells();
            for (NodeId n : {nodes...}) {
                assert(n.is_valid() && n.value < m_nodes.size() && "UnstructuredMesh::add_cell: node id out of range");
                m_c2n.push_back(n);
            }

            // 3. No neighbor is known yet across the new cell's faces; real adjacency will be
            // computed separately, see build_connectivity()
            for (SizeT i = 0; i < Cell::Faces; ++i) {
                m_c2c_raw.push_back(NeighborRef::NO_NEIGHBOR);
            }

            m_cell_variables.resize_all(nb_cells());
            return CellId{static_cast<CellId>(index)};
        }

        /**
         * @brief Computes cell-to-cell, face-to-face and face-to-cell adjacency from the mesh's
         * current set of nodes/edges/faces/cells.
         *
         * Must be called explicitly after the mesh has been fully populated (e.g. after reading a
         * file, or after manually calling add_edge()/add_face()/add_cell()); it is not run
         * automatically and is not incremental, so re-run it after further modifying the mesh.
         *
         * It never creates new Face/Edge entities: a geometric boundary facet of a cell without a
         * pre-existing Face for it (or a free edge of a face without a pre-existing Edge for it)
         * is left with a NeighborRef::Type::None neighbor rather than materializing one, matching
         * the convention that faces/edges only exist when explicitly added (e.g. tagged on read).
         *
         * @throw std::logic_error if the mesh is non-manifold: a facet shared by more than 2
         *        cells, or an edge shared by more than 2 faces.
         */
        void build_connectivity() {
            build_cell_connectivity();
            build_face_connectivity();
        }

        // -------------------------------------------------------------------------
        // Methods to access incident/adjacent cells
        // -------------------------------------------------------------------------
        /**
         * @brief Gets the position of a node.
         * @param id Id of the node to query.
         * @return Const reference to the node's position.
         * @pre @p id must be a valid id previously returned by add_node().
         */
        [[nodiscard]] const Point3d &node(NodeId id) const {
            assert(id.is_valid() && id.value < m_nodes.size() && "UnstructuredMesh::node: node id out of range");
            return m_nodes[id.value];
        }

        /**
         * @brief Gets a zero-copy view onto the node ids of a face.
         * @param id Id of the face to query.
         * @return A span of Face::Nodes node ids, in the order given to add_face().
         * @pre @p id must be a valid id previously returned by add_face().
         */
        [[nodiscard]] std::span<const NodeId> face_nodes(FaceId id) const {
            assert(id.is_valid() && id.value < nb_faces() && "UnstructuredMesh::face_nodes: face id out of range");
            return std::span<const NodeId>(m_f2n).subspan(static_cast<SizeT>(id.value) * Face::Nodes, Face::Nodes);
        }

        /**
         * @brief Gets a zero-copy view onto the node ids of a 3D cell.
         * @param id Id of the cell to query.
         * @return A span of Cell::Nodes node ids, in the order given to add_cell().
         * @pre @p id must be a valid id previously returned by add_cell().
         */
        [[nodiscard]] std::span<const NodeId> cell_nodes(CellId id) const {
            assert(id.is_valid() && id.value < nb_cells() && "UnstructuredMesh::cell_nodes: cell id out of range");
            return std::span<const NodeId>(m_c2n).subspan(static_cast<SizeT>(id.value) * Cell::Nodes, Cell::Nodes);
        }

        /**
         * @brief Gets a zero-copy view onto the two node ids of an edge.
         * @param id Id of the edge to query.
         * @return A span of the 2 node ids of the edge, in the order given to add_edge().
         * @pre @p id must be a valid id previously returned by add_edge().
         */
        [[nodiscard]] std::span<const NodeId> edge_nodes(EdgeId id) const {
            assert(id.is_valid() && id.value < nb_edges() && "UnstructuredMesh::edge_nodes: edge id out of range");
            return std::span<const NodeId>(m_e2n).subspan(static_cast<SizeT>(id.value) * 2, 2);
        }

        /**
         * @brief Gets the neighbor of a cell across one of its local faces.
         * @param id Id of the cell to query.
         * @param local_face Local face index in [0, Cell::Faces).
         * @return Cell3D (another cell, by index) if one exists across that face; otherwise
         *         Boundary2D (a Face, by index) if one was materialized for that boundary facet;
         *         otherwise None.
         * @pre @p id must be a valid id previously returned by add_cell().
         * @pre @p local_face must be in [0, Cell::Faces).
         * @pre build_connectivity() must have been called since the mesh was last modified.
         */
        [[nodiscard]] NeighborRef cell_neighbor(CellId id, SizeT local_face) const {
            assert(id.is_valid() && id.value < nb_cells() && "UnstructuredMesh::cell_neighbor: cell id out of range");
            assert(local_face < Cell::Faces && "UnstructuredMesh::cell_neighbor: local face index out of range");
            return NeighborRef::decode(m_c2c_raw[static_cast<SizeT>(id.value) * Cell::Faces + local_face]);
        }

        /**
         * @brief Gets the neighbor of a face across one of its local edges.
         * @param id Id of the face to query.
         * @param local_edge Local edge index in [0, Face::Edges).
         * @return Cell3D (another face, by index) if one exists across that edge; otherwise
         *         Boundary2D (an Edge, by index) if one was materialized for that free edge;
         *         otherwise None.
         * @pre @p id must be a valid id previously returned by add_face().
         * @pre @p local_edge must be in [0, Face::Edges).
         * @pre build_connectivity() must have been called since the mesh was last modified.
         */
        [[nodiscard]] NeighborRef face_neighbor(FaceId id, SizeT local_edge) const {
            assert(id.is_valid() && id.value < nb_faces() && "UnstructuredMesh::face_neighbor: face id out of range");
            assert(local_edge < Face::Edges && "UnstructuredMesh::face_neighbor: local edge index out of range");
            return NeighborRef::decode(m_f2f_raw[static_cast<SizeT>(id.value) * Face::Edges + local_edge]);
        }

        /**
         * @brief Gets the cell(s) adjacent to a face.
         * @param id Id of the face to query.
         * @return The 1 or 2 cells bordering the face, in no particular order; an unused slot (a
         *         boundary face borders only 1 cell) holds CellId::invalid_id.
         * @pre @p id must be a valid id previously returned by add_face().
         * @pre build_connectivity() must have been called since the mesh was last modified.
         */
        [[nodiscard]] std::array<CellId, 2> face_cells(FaceId id) const {
            assert(id.is_valid() && id.value < nb_faces() && "UnstructuredMesh::face_cells: face id out of range");
            const auto base = static_cast<SizeT>(id.value) * 2;
            return {m_f2c_raw[base], m_f2c_raw[base + 1]};
        }

        // -------------------------------------------------------------------------
        // Methods to handle variables on the mesh cells
        // -------------------------------------------------------------------------
        /**
         * @brief Adds a variable of type @p TVar and name @p name on the cells of type @p TCell
         * (node, edge, face or 3D cell).
         * @tparam TVar Type of the variable to create.
         * @tparam TCell Cell type the variable is defined on.
         * @param name Name of the variable.
         * @return A reference onto the newly created variable.
         * @pre TCell must be one of CellType::Node, CellType::Edge, CellType::Face or CellType::Cell.
         * @throw std::invalid_argument if a variable of that name already exists.
         */
        template<typename TVar, CellType TCell>
        Variable<TVar> &add_variable(const std::string &name) {
            static_assert(TCell == CellType::Node || TCell == CellType::Edge || TCell == CellType::Face ||
                              TCell == CellType::Cell,
                          "add_variable() : TCell must be Node, Edge, Face or Cell !");
            if constexpr (TCell == CellType::Node) {
                return m_node_variables.add<TVar>(name);
            } else if constexpr (TCell == CellType::Edge) {
                return m_edge_variables.add<TVar>(name);
            } else if constexpr (TCell == CellType::Face) {
                return m_face_variables.add<TVar>(name);
            } else if constexpr (TCell == CellType::Cell) {
                return m_cell_variables.add<TVar>(name);
            }
        }

        /**
         * @brief Removes the variable of name @p name attached to the cells of type @p TCell.
         * @tparam TCell Cell type the variable is defined on.
         * @param name Name of the variable to remove.
         * @return true if the variable was found and removed, false otherwise.
         * @pre TCell must be one of CellType::Node, CellType::Edge, CellType::Face or CellType::Cell.
         */
        template<CellType TCell>
        bool remove_variable(const std::string &name) {
            static_assert(TCell == CellType::Node || TCell == CellType::Edge || TCell == CellType::Face ||
                              TCell == CellType::Cell,
                          "remove_variable() : TCell must be Node, Edge, Face or Cell !");

            if constexpr (TCell == CellType::Node) {
                return m_node_variables.remove(name);
            } else if constexpr (TCell == CellType::Edge) {
                return m_edge_variables.remove(name);
            } else if constexpr (TCell == CellType::Face) {
                return m_face_variables.remove(name);
            } else if constexpr (TCell == CellType::Cell) {
                return m_cell_variables.remove(name);
            }
        }

        /**
         * @brief Checks whether a variable of name @p AName exists on the cells of type @tparam TCell.
         * @tparam TCell Cell type the variable would be defined on.
         * @param AName Name of the variable to look for.
         * @return true if the variable exists, false otherwise.
         * @pre TCell must be one of CellType::Node, CellType::Edge, CellType::Face or CellType::Cell.
         */
        template<CellType TCell>
        bool has_variable(const std::string &AName) {
            static_assert(TCell == CellType::Node || TCell == CellType::Edge || TCell == CellType::Face ||
                              TCell == CellType::Cell,
                          "has_variable() : TCell must be Node, Edge, Face or Cell !");

            if constexpr (TCell == CellType::Node) {
                return m_node_variables.has(AName);
            } else if constexpr (TCell == CellType::Edge) {
                return m_edge_variables.has(AName);
            } else if constexpr (TCell == CellType::Face) {
                return m_face_variables.has(AName);
            } else if constexpr (TCell == CellType::Cell) {
                return m_cell_variables.has(AName);
            }
        }

        /** @copydoc has_variable(const std::string&) */
        template<CellType TCell>
        bool has_variable(const std::string &AName) const {
            static_assert(TCell == CellType::Node || TCell == CellType::Edge || TCell == CellType::Face ||
                              TCell == CellType::Cell,
                          "has_variable() : TCell must be Node, Edge, Face or Cell !");

            if constexpr (TCell == CellType::Node) {
                return m_node_variables.has(AName);
            } else if constexpr (TCell == CellType::Edge) {
                return m_edge_variables.has(AName);
            } else if constexpr (TCell == CellType::Face) {
                return m_face_variables.has(AName);
            } else if constexpr (TCell == CellType::Cell) {
                return m_cell_variables.has(AName);
            }
        }

        /**
         * @brief Gets the variable of name @p name.
         * @tparam TVar Type of the variable's items.
         * @tparam TCell Cell type the variable is defined on.
         * @param name Name of the variable to retrieve.
         * @return A reference onto the variable.
         * @pre TCell must be one of CellType::Node, CellType::Edge, CellType::Face or CellType::Cell.
         * @throw when no variable of name @p name exists for the given cell type.
         */
        template<typename TVar, CellType TCell>
        Variable<TVar> &get_variable(const std::string &name) {
            static_assert(TCell == CellType::Node || TCell == CellType::Edge || TCell == CellType::Face ||
                              TCell == CellType::Cell,
                          "get_variable() : TCell must be Node, Edge, Face or Cell !");

            if constexpr (TCell == CellType::Node) {
                return m_node_variables.get<TVar>(name);
            } else if constexpr (TCell == CellType::Edge) {
                return m_edge_variables.get<TVar>(name);
            } else if constexpr (TCell == CellType::Face) {
                return m_face_variables.get<TVar>(name);
            } else if constexpr (TCell == CellType::Cell) {
                return m_cell_variables.get<TVar>(name);
            }
        }

        /** @copydoc get_variable(const std::string&) */
        template<typename TVar, CellType TCell>
        const Variable<TVar> &get_variable(const std::string &name) const {
            static_assert(TCell == CellType::Node || TCell == CellType::Edge || TCell == CellType::Face ||
                              TCell == CellType::Cell,
                          "get_variable() : TCell must be Node, Edge, Face or Cell !");

            if constexpr (TCell == CellType::Node) {
                return m_node_variables.get<TVar>(name);
            } else if constexpr (TCell == CellType::Edge) {
                return m_edge_variables.get<TVar>(name);
            } else if constexpr (TCell == CellType::Face) {
                return m_face_variables.get<TVar>(name);
            } else if constexpr (TCell == CellType::Cell) {
                return m_cell_variables.get<TVar>(name);
            }
        }

    private:
        /**
         * @brief Builds a sorted key of node values from a facet's local node indices.
         * @tparam N Number of local node indices (i.e. nodes of the facet).
         * @param verts The parent entity's node ids (e.g. cell_nodes()/face_nodes()).
         * @param local Local node indices (into @p verts) forming the facet.
         * @return The sorted node values of the facet, used as a facet identity key.
         */
        template<std::size_t N>
        static std::vector<UInt> facet_key(std::span<const NodeId> verts, const std::array<SizeT, N> &local) {
            std::vector<UInt> key;
            key.reserve(N);
            for (SizeT l : local) {
                key.push_back(verts[l].value);
            }
            std::sort(key.begin(), key.end());
            return key;
        }

        /**
         * @brief Computes cell-to-cell adjacency (m_c2c_raw) and the cell side of face-to-cell
         * adjacency (m_f2c_raw), see build_connectivity().
         * @throw std::logic_error if a facet is shared by more than 2 cells.
         */
        void build_cell_connectivity() {
            // Existing faces, indexed by their node-set, so a matched boundary/tagged-internal
            // facet can be linked to the Face materialized for it (none is ever created here).
            std::map<std::vector<UInt>, FaceId> face_by_key;
            for (UInt f = 0; f < nb_faces(); ++f) {
                auto verts = face_nodes(FaceId{f});
                std::vector<UInt> key;
                key.reserve(verts.size());
                for (auto v : verts) {
                    key.push_back(v.value);
                }
                std::sort(key.begin(), key.end());
                face_by_key[key] = FaceId{f};
            }

            std::fill(m_f2c_raw.begin(), m_f2c_raw.end(), CellId::invalid_id);

            // Every (cell, local face) touching each facet key.
            std::map<std::vector<UInt>, std::vector<std::pair<CellId, SizeT>>> touches;
            for (UInt c = 0; c < nb_cells(); ++c) {
                auto verts = cell_nodes(CellId{c});
                for (SizeT lf = 0; lf < Cell::Faces; ++lf) {
                    touches[facet_key(verts, Cell::FaceNodes[lf])].emplace_back(CellId{c}, lf);
                }
            }

            for (auto &[key, list] : touches) {
                if (list.size() == 1) {
                    const auto [cell, lf] = list[0];
                    if (auto it = face_by_key.find(key); it != face_by_key.end()) {
                        m_c2c_raw[cell.value * Cell::Faces + lf] =
                            NeighborRef::encode_boundary(static_cast<Int>(it->second.value));
                        m_f2c_raw[it->second.value * 2] = cell;
                    } else {
                        m_c2c_raw[cell.value * Cell::Faces + lf] = NeighborRef::NO_NEIGHBOR;
                    }
                } else if (list.size() == 2) {
                    const auto [cellA, lfA] = list[0];
                    const auto [cellB, lfB] = list[1];
                    m_c2c_raw[cellA.value * Cell::Faces + lfA] = static_cast<Int>(cellB.value);
                    m_c2c_raw[cellB.value * Cell::Faces + lfB] = static_cast<Int>(cellA.value);
                    if (auto it = face_by_key.find(key); it != face_by_key.end()) {
                        m_f2c_raw[it->second.value * 2] = cellA;
                        m_f2c_raw[it->second.value * 2 + 1] = cellB;
                    }
                } else {
                    throw std::logic_error(
                        "UnstructuredMesh::build_connectivity: non-manifold mesh (a facet is shared by more "
                        "than 2 cells)");
                }
            }
        }

        /**
         * @brief Computes face-to-face adjacency (m_f2f_raw), see build_connectivity().
         * @throw std::logic_error if an edge is shared by more than 2 faces.
         */
        void build_face_connectivity() {
            // Existing edges, indexed by their node-set, so a matched free edge can be linked to
            // the Edge materialized for it (none is ever created here).
            std::map<std::vector<UInt>, EdgeId> edge_by_key;
            for (UInt e = 0; e < nb_edges(); ++e) {
                auto n = edge_nodes(EdgeId{e});
                std::vector<UInt> key{n[0].value, n[1].value};
                std::sort(key.begin(), key.end());
                edge_by_key[key] = EdgeId{e};
            }

            // Every (face, local edge) touching each facet key.
            std::map<std::vector<UInt>, std::vector<std::pair<FaceId, SizeT>>> touches;
            for (UInt f = 0; f < nb_faces(); ++f) {
                auto verts = face_nodes(FaceId{f});
                for (SizeT le = 0; le < Face::Edges; ++le) {
                    touches[facet_key(verts, Face::EdgeNodes[le])].emplace_back(FaceId{f}, le);
                }
            }

            for (auto &[key, list] : touches) {
                if (list.size() == 1) {
                    const auto [face, le] = list[0];
                    if (auto it = edge_by_key.find(key); it != edge_by_key.end()) {
                        m_f2f_raw[face.value * Face::Edges + le] =
                            NeighborRef::encode_boundary(static_cast<Int>(it->second.value));
                    } else {
                        m_f2f_raw[face.value * Face::Edges + le] = NeighborRef::NO_NEIGHBOR;
                    }
                } else if (list.size() == 2) {
                    const auto [faceA, leA] = list[0];
                    const auto [faceB, leB] = list[1];
                    m_f2f_raw[faceA.value * Face::Edges + leA] = static_cast<Int>(faceB.value);
                    m_f2f_raw[faceB.value * Face::Edges + leB] = static_cast<Int>(faceA.value);
                } else {
                    throw std::logic_error(
                        "UnstructuredMesh::build_connectivity: non-manifold face set (an edge is shared by "
                        "more than 2 faces)");
                }
            }
        }

        // nodes location
        std::vector<Point3d> m_nodes;

        // --- 3D Cells (Tet / Hex) ---
        std::vector<NodeId> m_c2n; // N_3d * Cell::Nodes
        // N_3d * Cell::Faces, signed neighbor code per local face: >= 0 is another cell's index,
        // < 0 (and != NeighborRef::NO_NEIGHBOR) is an encoded boundary FaceId,
        // NeighborRef::NO_NEIGHBOR is no neighbor at all (see NeighborRef / build_connectivity()).
        std::vector<Int> m_c2c_raw;

        // --- 2D/3D Faces (Tri / Quad) ---
        std::vector<NodeId> m_f2n; // N_2d * Face::Nodes
        // N_2d * Face::Edges, signed neighbor code per local edge, same convention as m_c2c_raw
        // one dimension down (>= 0: another face, encoded: a boundary EdgeId,
        // NeighborRef::NO_NEIGHBOR: none).
        std::vector<Int> m_f2f_raw;
        std::vector<CellId> m_f2c_raw; // N_2d * 2 (the 1 or 2 adjacent 3D cells; CellId::invalid_id if none)

        // --- Ridge cells (Edges) ---
        std::vector<NodeId> m_e2n;

        // Variable registries, one per type of cells
        VariableRegistry m_node_variables;
        VariableRegistry m_edge_variables;
        VariableRegistry m_face_variables;
        VariableRegistry m_cell_variables;
    };

    // -----------------------------------------------------------------------------
    using SimplicialMesh = UnstructuredMesh<SimplicialTraits>; // Triangles + Tets
    using CubicMesh = UnstructuredMesh<CubicTraits>;           // Quads + Hexes
    // -----------------------------------------------------------------------------
} // namespace gecko
