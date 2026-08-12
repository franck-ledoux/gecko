/*----------------------------------------------------------------------------*/
#pragma once
#include <vector>
//#include <mdspan>
#include <cstdint>
#include <expected>
// -----------------------------------------------------------------------------
#include <gecko/math/Point3d.h>
#include <gecko/mesh/VariableRegistry.h>
#include <gecko/utils/Types.h>
// -----------------------------------------------------------------------------
//namespace stdex = std::experimental;
// -----------------------------------------------------------------------------
// Traits to define mesh topology
// -----------------------------------------------------------------------------
namespace gecko {
    /** @brief Topology traits describing simplicial cells (tetrahedra / triangles). */
    struct SimplicialTraits {
        /** @brief Topology of a 3D cell (tetrahedron): 4 nodes, 4 faces. */
        struct Cell {
            static constexpr size_t Nodes = 4; ///< Number of nodes per cell.
            static constexpr size_t Faces = 4; ///< Number of faces per cell.
        };

        /** @brief Topology of a 2D face (triangle): 3 nodes, 3 edges. */
        struct Face {
            static constexpr size_t Nodes = 3; ///< Number of nodes per face.
            static constexpr size_t Edges = 3; ///< Number of edges per face.
        };
    };

    // -----------------------------------------------------------------------------
    /** @brief Topology traits describing cubic cells (hexahedra / quadrilaterals). */
    struct CubicTraits {
        /** @brief Topology of a 3D cell (hexahedron): 8 nodes, 6 faces. */
        struct Cell {
            static constexpr size_t Nodes = 8; ///< Number of nodes per cell.
            static constexpr size_t Faces = 6; ///< Number of faces per cell.
        };

        /** @brief Topology of a 2D face (quadrilateral): 4 nodes, 4 edges. */
        struct Face {
            static constexpr size_t Nodes = 4; ///< Number of nodes per face.
            static constexpr size_t Edges = 4; ///< Number of edges per face.
        };
    };

    /*----------------------------------------------------------------------------*/
    //template <typename T, size_t Extent>
    //using View2D = stdex::mdspan<T, stdex::extents<size_t, stdex::dynamic_extent, Extent>>;
    /*----------------------------------------------------------------------------*/
    /** @brief Dimension of a mesh cell type. */
    enum class CellType { Node, Edge, Face, Cell };
    /** @brief Placeholder for future 3D cell data (currently unused). */
    struct CellData {};

    /*----------------------------------------------------------------------------*/
    /** @brief Encoded reference to a mesh cell's neighbor, either a 3D cell or a 2D boundary face. */
    struct NeighborRef {
        /** @brief Kind of entity a NeighborRef points to. */
        enum class Type { Cell3D, Boundary2D };

        Type type;     ///< Whether the referenced entity is a 3D cell or a 2D boundary face.
        int32_t index; ///< Index of the referenced entity within its own container.
    };

    /**
     * @brief Encodes a boundary face index into the signed neighbor code convention.
     * @param boundary_idx Index of the boundary face (>= 0).
     * @return A negative code identifying the boundary face, decodable by decode_neighbor().
     */
    /*----------------------------------------------------------------------------*/
    inline int32_t encode_boundary(int32_t boundary_idx) { return -boundary_idx - 1; }
    /*----------------------------------------------------------------------------*/

    /**
     * @brief Decodes a signed neighbor code into a NeighborRef.
     *
     * A non-negative code refers to a 3D cell index; a negative code refers to
     * an encoded boundary face index (see encode_boundary()).
     *
     * @param code Signed neighbor code to decode.
     * @return The decoded NeighborRef, tagged as Cell3D or Boundary2D.
     */
    inline NeighborRef decode_neighbor(int32_t code) {
        if (code >= 0) {
            return {NeighborRef::Type::Cell3D, code};
        } else {
            return {NeighborRef::Type::Boundary2D, -code - 1};
        }
    }

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
        explicit UnstructuredMesh(size_t capacity = 64) { m_nodes.reserve(capacity); }

        // Vues multidimensionnelles (Zero-overhead C++23)
        /*  auto cell_nodes()  { return View2D<int32_t, Cell::Vertices>(m_c2n_raw.data(), m_c2n_raw.size() / Cell::Nodes); }
          auto cell_neighbors() { return View2D<int32_t, Cell::Faces>(m_c2c_raw.data(), m_c2c_raw.size() / Cell::Faces); }

          auto face_nodes()  { return View2D<int32_t, Face::Vertices>(m_f2n_raw.data(), m_f2n_raw.size() / Face::Nodes); }
          auto face_neighbors() { return View2D<int32_t, Face::Edges>(m_f2f_raw.data(), m_f2f_raw.size() / Face::Edges); }
          auto face_cells()   { return View2D<int32_t, 2>(m_f2c_raw.data(), m_f2c_raw.size() / 2); }
      */

        /** @brief Gets the number of nodes stored in the mesh. @return The node count. */
        [[nodiscard]] UInt nb_nodes() const { return m_nodes.size(); }
        /** @brief Gets the number of edges stored in the mesh. @return The edge count. */
        [[nodiscard]] UInt nb_edges() const { return m_e2n.size() / 2; }
        /** @brief Gets the number of faces stored in the mesh. @return The face count. */
        [[nodiscard]] UInt nb_faces() const { return m_c2n.size() / Cell::Vertices; }
        /** @brief Gets the number of 3D cells stored in the mesh. @return The cell count. */
        [[nodiscard]] UInt nb_cells() const { return m_f2n.size() / Face::Vertices; }

        /**
         * @brief Adds a new node to the mesh at the given position.
         * @param pnt Position of the new node.
         * @return The id of the newly created node.
         */
        [[nodiscard]] NodeId add_node(const Point3d &pnt) {
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
        // Surcharge pratique par coordonnées directes (ex: x, y, z)
        [[nodiscard]] NodeId add_node(Float x, Float y, Float z = 0.0) { return add_node(Point3d{x, y, z}); }

        /*    template <typename... Args>
            size_t new_face(Args... ANodes) {
                // 1. Contrôle à la compilation
                static_assert(sizeof...(Args) == Face::Nodes,
                    "Erreur : Le nombre de nœuds ne correspond pas à la topologie de la face !");

                // 2. Contrôle optionnel des types (s'assurer qu'on ne passe que des entiers/indices)
                static_assert((std::is_convertible_v<Args, size_t> && ...),
                    "Erreur : Tous les arguments doivent être convertibles en size_t !");

                // 3. Construction directe de la face sans std::copy
                FaceType face{ static_cast<size_t>(nodes)... };

                m_faces.push_back(face);
                m_face_registry.resize_all(m_faces.size());
                return m_faces.size() - 1;
            }*/
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

        /**
         * @brief Gets the underlying `std::vector<TVar>` storing a variable's values.
         * @tparam TVar Type of the variable's items.
         * @tparam TCell Cell type the variable is defined on.
         * @param name Name of the variable to retrieve.
         * @return A reference onto the variable's data vector.
         * @pre TCell must be one of CellType::Node, CellType::Edge, CellType::Face or CellType::Cell.
         * @throw when no variable of name @p name exists for the given cell type.
         */
        template<typename TVar, CellType TCell>
        std::vector<TVar> &get_variable(const std::string &name) {
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
        // nodes location
        std::vector<Point3d> m_nodes;

        // --- 3D Cells (Tet / Hex) ---
        std::vector<Int> m_c2n;     // N_3d * Cell::Nodes
        std::vector<Int> m_c2c_raw; // N_3d * Cell::Faces (signé : + si 3D, - si 2D)

        // --- Élément Surface (Tri / Quad) ---
        std::vector<Int> m_f2n;     // N_2d * Face::Nodes
        std::vector<Int> m_f2f_raw; // N_2d * Face::Edges (voisins 2D par arête)
        std::vector<Int> m_f2c_raw; // N_2d * 2 (indices des 2 éléments 3D adjacents)

        // --- Élément Ridge (Edge) ---
        std::vector<Int> m_e2n;

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
