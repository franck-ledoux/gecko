#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#include <gecko/block/Blocking.h>
#include <gecko/math/BezierCurve.h>
#include <gecko/math/Point3d.h>

#include "GeomModelFacade.h"

namespace gecko::python {

    /**
     * @class BlockingFacade
     * @brief Python-facing façade over gecko::Blocking<FacetedGeometry, TEdgeCurve>: every method
     * takes/returns only primitive types — corners as lists of (x,y,z) tuples, faces/blocks as
     * plain int ids assigned by this façade, never a gecko::Point3d or CGAL dart/attribute handle
     * — see docs/user-guide/python.md.
     *
     * TEdgeCurve's degree (1 = straight, 3 = cubic Bezier) is chosen at construction time via a
     * plain int, instead of being a compile-time Python-visible axis: internally this holds one of
     * a small std::variant of degree-specific implementations, matching the 2 instantiations used
     * across block/tst/ (see e.g. blocking_curve_evaluation_tests.cpp).
     */
    class BlockingFacade {
    public:
        /**
         * @brief Constructor.
         * @param model The geometric model to build against; must outlive this Blocking (enforced
         *        on the Python side via py::keep_alive).
         * @param degree Edge curve degree: 1 for straight edges, 3 for cubic Bezier edges.
         * @throw std::invalid_argument if @p degree is neither 1 nor 3.
         */
        explicit BlockingFacade(const GeomModelFacade &model, int degree = 1);

        /**
         * @brief Creates a new, unsewn quad block.
         * @param corners The 4 corner positions, each an (x,y,z) triple, in perimeter order.
         * @return The id assigned to the newly created face.
         * @throw std::invalid_argument if @p corners doesn't have exactly 4 entries.
         */
        int create_quad_block(const std::vector<std::array<double, 3>> &corners);

        /**
         * @brief Creates a new, unsewn hex block.
         * @param corners The 8 corner positions, each an (x,y,z) triple — see
         *        gecko::Blocking::create_hex_block for the expected HEX8 ordering.
         * @return The id assigned to the newly created block.
         * @throw std::invalid_argument if @p corners doesn't have exactly 8 entries.
         */
        int create_hex_block(const std::vector<std::array<double, 3>> &corners);

        /** @brief Auto-detects and sews coincident blocks created so far. Not incremental. */
        void build_connectivity();

        /**
         * @brief Classifies every cell onto the geometric model and refits geometry accordingly.
         * @param tol_vertex Tolerance for snapping onto a vertex.
         * @param tol_curve_surface Tolerance for snapping onto a curve/surface/volume; defaults to
         *        @p tol_vertex when left at -1.0.
         */
        void classify(double tol_vertex, double tol_curve_surface = -1.0);

        /**
         * @brief Gets the number of cells of a given dimension in the blocking.
         * @param dim Cell dimension, in [0,3].
         * @return The cell count.
         * @throw std::invalid_argument if @p dim is not in [0,3].
         */
        [[nodiscard]] std::size_t nb_cells(int dim) const;

        /** @brief Checks the topological validity of the underlying map. @return true if valid. */
        [[nodiscard]] bool is_valid_topology() const;
        /** @brief Checks whether this blocking has no 3-cell (hex block). @return true if purely 2D. */
        [[nodiscard]] bool is_purely_2d() const;
        /**
         * @brief Checks whether a face can be deleted.
         * @param face_id A face id returned by create_quad_block() (or a hex block's boundary
         *        face — not currently reachable through this façade).
         * @return true if the face can be deleted.
         * @throw std::out_of_range if @p face_id is not a known face id.
         */
        [[nodiscard]] bool can_delete_face(int face_id) const;
        /**
         * @brief Deletes a face from the structure.
         * @param face_id A face id returned by create_quad_block().
         * @throw std::out_of_range if @p face_id is not a known face id.
         */
        void delete_face(int face_id);

        /**
         * @brief Generates a quad/hex mesh reproducing the blocking's geometry and writes it to a
         * VTK legacy ASCII file.
         * @param subdivisions Number of intervals to subdivide every parametric axis into (>= 1).
         * @param path Output file path.
         */
        void write_vtk(int subdivisions, const std::string &path);

        /** @brief Gets the id of every corner node of the block structure. @return The node ids. */
        [[nodiscard]] std::vector<int> node_ids() const;
        /**
         * @brief Gets the position of a corner node.
         * @param node_id A node id, from node_ids().
         * @return Its (x,y,z) position.
         * @throw std::out_of_range if @p node_id is not a known node id.
         */
        [[nodiscard]] std::array<double, 3> node_position(int node_id) const;
        /**
         * @brief Moves a corner node, refitting every incident edge/face/block geometry.
         * @param node_id A node id, from node_ids().
         * @param x New X coordinate.
         * @param y New Y coordinate.
         * @param z New Z coordinate.
         * @throw std::out_of_range if @p node_id is not a known node id.
         */
        void move_node(int node_id, double x, double y, double z);

        /**
         * @brief Generates the mesh reproducing the blocking and returns its node positions.
         * @param subdivisions Number of intervals to subdivide every parametric axis into (>= 1).
         * @return One (x,y,z) triple per mesh node.
         */
        [[nodiscard]] std::vector<std::array<double, 3>> mesh_vertices(int subdivisions);
        /**
         * @brief Generates the mesh reproducing the blocking and returns its quads.
         * @param subdivisions Number of intervals to subdivide every parametric axis into (>= 1).
         * @return One quadruple of mesh_vertices() indices per quad.
         */
        [[nodiscard]] std::vector<std::array<int, 4>> mesh_quads(int subdivisions);
        /**
         * @brief Generates the mesh reproducing the blocking and returns its hexes.
         * @param subdivisions Number of intervals to subdivide every parametric axis into (>= 1).
         * @return One 8-tuple of mesh_vertices() indices per hex.
         */
        [[nodiscard]] std::vector<std::array<int, 8>> mesh_hexes(int subdivisions);

        /** @brief Per-degree state; public only so the .cpp's free-function helpers can name it —
         * never part of the class' actual (Python-facing) interface. */
        template<std::size_t N>
        struct Impl {
            using BlockingT = Blocking<FacetedGeometry, BezierCurve<N, Point3d>>;

            BlockingT blocking;
            std::unordered_map<int, typename BlockingT::Face> faces_by_id;
            std::unordered_map<int, typename BlockingT::Node> nodes_by_id;
            int next_face_id = 0;
            int next_block_id = 0;
            int next_node_id = 0;

            explicit Impl(const FacetedGeometry &geom) : blocking(geom) {}

            /** @brief Registers every not-yet-known corner node of the structure, so a freshly
             * created block's corners become addressable by id. Cheap and idempotent: nodes already
             * present keep their id. */
            void index_new_nodes() {
                auto &map = blocking.cmap();
                for (auto it = map.template attributes<0>().begin(), itend = map.template attributes<0>().end();
                     it != itend;
                     ++it) {
                    const bool known = std::any_of(nodes_by_id.begin(), nodes_by_id.end(), [&it](const auto &entry) {
                        return entry.second == it;
                    });
                    if (!known) nodes_by_id.emplace(next_node_id++, it);
                }
            }
        };

    private:
        std::variant<Impl<1>, Impl<3>> m_impl;
    };

} // namespace gecko::python
