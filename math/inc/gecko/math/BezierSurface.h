#pragma once

#include <array>
#include <cassert>
#include <cstddef>

namespace gecko {
    /**
     * @class BezierSurface
     * @brief Represents a tensor-product Bezier surface of degree @p TN x @p TN, with a
     * `(TN+1) x (TN+1)` grid of control points of type @p TPointT.
     *
     * Generalizes BezierCurve's De Casteljau evaluation to 2 parametric dimensions: evaluating at
     * `(u,v)` reduces each of the `TN+1` control-point rows along `v`, then reduces the resulting
     * `TN+1` points along `u`.
     *
     * @tparam TN The degree of the surface in each parametric direction (control grid is
     * `(TN+1) x (TN+1)`).
     * @tparam TPointT The geometric point type, following gecko's point/vector affine-space split
     * (e.g. gecko::Point3d): `TPointT - TPointT` yields a difference/vector type `V`, `V * double`
     * yields `V`, and `TPointT + V` yields `TPointT` (see BezierCurve.h).
     */
    template<std::size_t TN, typename TPointT>
    class BezierSurface {
    public:
        /** @brief The degree of the surface in each parametric direction. */
        static constexpr std::size_t Degree = TN;

        /** @brief The number of control points along each parametric direction (N + 1). */
        static constexpr std::size_t NumControlPoints = TN + 1;

        /** @brief Row-major control point grid: `m_control_points[i][j]`, `i` along `u`, `j` along `v`. */
        using ControlGrid = std::array<std::array<TPointT, NumControlPoints>, NumControlPoints>;

        /** @brief Default constructor. Every control point default-initializes to `TPointT{}`. */
        BezierSurface() = default;

        /**
         * @brief Constructor. Constructs a Bezier surface from an explicit control grid.
         * @param AControlPoints A `(TN+1) x (TN+1)` grid of control points, `[i][j]` along `u`/`v`.
         */
        explicit BezierSurface(const ControlGrid &AControlPoints) : m_control_points(AControlPoints) {}

        /**
         * @brief Gets a read-only reference to the control grid.
         * @return Const reference to the `(TN+1) x (TN+1)` control point grid.
         */
        const ControlGrid &control_points() const { return m_control_points; }

        /**
         * @brief Accesses a control point by its 2D grid index.
         * @param AI Index along `u` (0 to N).
         * @param AJ Index along `v` (0 to N).
         * @return Reference to the requested control point.
         * @note Named `control_point`, not `operator()`, to avoid an unresolvable overload
         *       ambiguity with the `(double,double)` evaluation `operator()` below for calls with
         *       integer-literal arguments (equally convertible to `std::size_t` or `double`).
         */
        TPointT &control_point(std::size_t AI, std::size_t AJ) { return m_control_points[AI][AJ]; }

        /**
         * @brief Accesses a control point by its 2D grid index.
         * @param AI Index along `u` (0 to N).
         * @param AJ Index along `v` (0 to N).
         * @return Const reference to the requested control point.
         */
        const TPointT &control_point(std::size_t AI, std::size_t AJ) const { return m_control_points[AI][AJ]; }

        /**
         * @brief Evaluates the point on the surface at parameters (u, v) using De Casteljau's
         * algorithm applied first along `v`, then along `u`.
         * @param AU Parametric coordinate along `u`, in the range [0.0, 1.0].
         * @param AV Parametric coordinate along `v`, in the range [0.0, 1.0].
         * @return The evaluated point on the surface.
         * @pre Both AU and AV must satisfy 0.0 <= x <= 1.0 (enforced via assert in Debug builds).
         */
        TPointT value(double AU, double AV) const {
            assert(AU >= 0.0 && AU <= 1.0 && "BezierSurface::value: AU must be in the range [0.0, 1.0]");
            assert(AV >= 0.0 && AV <= 1.0 && "BezierSurface::value: AV must be in the range [0.0, 1.0]");
            std::array<TPointT, NumControlPoints> reduced_along_v;
            for (std::size_t i = 0; i < NumControlPoints; ++i) {
                reduced_along_v[i] = de_casteljau(m_control_points[i], AV);
            }
            return de_casteljau(reduced_along_v, AU);
        }

        /**
         * @brief Evaluates the surface at parameters (u, v) (functor operator shortcut).
         * @param AU Parametric coordinate along `u`.
         * @param AV Parametric coordinate along `v`.
         * @return The evaluated point on the surface.
         * @see value()
         */
        TPointT operator()(double AU, double AV) const { return value(AU, AV); }

    private:
        /** @brief The `(TN+1) x (TN+1)` grid storing the control points. */
        ControlGrid m_control_points{};

        /**
         * @brief Reduces a row of `NumControlPoints` points to a single point via 1D De Casteljau.
         * @param APoints The row of control points to reduce (passed by value: worked on in place).
         * @param AT Interpolation parameter in [0.0, 1.0].
         * @return The reduced point.
         */
        static TPointT de_casteljau(std::array<TPointT, NumControlPoints> APoints, double AT) {
            for (std::size_t k = 1; k <= TN; ++k) {
                for (std::size_t i = 0; i <= TN - k; ++i) {
                    APoints[i] = APoints[i] + (APoints[i + 1] - APoints[i]) * AT;
                }
            }
            return APoints[0];
        }
    };
} // namespace gecko
