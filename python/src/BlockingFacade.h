#pragma once

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

        /** @brief Per-degree state; public only so the .cpp's free-function helpers can name it —
         * never part of the class' actual (Python-facing) interface. */
        template<std::size_t N>
        struct Impl {
            Blocking<FacetedGeometry, BezierCurve<N, Point3d>> blocking;
            std::unordered_map<int, typename Blocking<FacetedGeometry, BezierCurve<N, Point3d>>::Face> faces_by_id;
            int next_face_id = 0;
            int next_block_id = 0;

            explicit Impl(const FacetedGeometry &geom) : blocking(geom) {}
        };

    private:
        std::variant<Impl<1>, Impl<3>> m_impl;
    };

} // namespace gecko::python
