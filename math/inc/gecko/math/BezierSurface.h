#pragma once

#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

#include <gecko/math/BezierCurve.h>

namespace gecko {
    /**
     * @class BezierSurface
     * @brief A tensor-product Bezier surface of any degree, with a `(degree+1) x (degree+1)` grid of
     * control points of type @p TPointT.
     *
     * The degree is a property of the object rather than of its type, for the reason given on
     * `BezierCurve`: a block structure has to be able to change order in place.
     *
     * @tparam TPointT The geometric point type, following gecko's point/vector affine-space split
     *                (e.g. gecko::Point3d).
     */
    template<typename TPointT>
    class BezierSurface {
    public:
        /** @brief Row-major control point grid: `grid[i][j]`, `i` along `u`, `j` along `v`. */
        using ControlGrid = std::vector<std::vector<TPointT>>;

        /**
         * @brief Builds an empty control grid for a given degree, every point default-initialized.
         * @param ADegree The degree in each parametric direction.
         * @return A `(ADegree+1) x (ADegree+1)` grid.
         */
        static ControlGrid make_grid(std::size_t ADegree) {
            return ControlGrid(ADegree + 1, std::vector<TPointT>(ADegree + 1));
        }

        /** @brief Default constructor: a degree-1 surface with every control point default-initialized. */
        BezierSurface() : m_control_points(make_grid(1)) {}

        /**
         * @brief Constructs a surface of a given degree, control points default-initialized.
         * @param ADegree The degree in each parametric direction.
         */
        explicit BezierSurface(std::size_t ADegree) : m_control_points(make_grid(ADegree)) {}

        /**
         * @brief Constructs a Bezier surface from an explicit control grid; the degree follows from
         * its size.
         * @param AControlPoints A square grid of control points, `[i][j]` along `u`/`v`.
         */
        explicit BezierSurface(ControlGrid AControlPoints) : m_control_points(std::move(AControlPoints)) {
            assert(m_control_points.size() >= 2 && m_control_points.size() == m_control_points.front().size() &&
                   "BezierSurface: the control grid must be square and at least 2x2");
        }

        /** @brief Gets the surface's degree in each parametric direction. @return The degree. */
        [[nodiscard]] std::size_t degree() const { return m_control_points.size() - 1; }

        /** @brief Gets the control point count along each direction. @return `degree() + 1`. */
        [[nodiscard]] std::size_t nb_control_points() const { return m_control_points.size(); }

        /**
         * @brief Gets a read-only reference to the control grid.
         * @return Const reference to the control point grid.
         */
        const ControlGrid &control_points() const { return m_control_points; }

        /**
         * @brief Accesses one control point.
         * @param AI Index along `u`, in `[0, degree()]`.
         * @param AJ Index along `v`, in `[0, degree()]`.
         * @return Reference to that control point.
         */
        TPointT &control_point(std::size_t AI, std::size_t AJ) { return m_control_points[AI][AJ]; }

        /**
         * @brief Accesses one control point.
         * @param AI Index along `u`, in `[0, degree()]`.
         * @param AJ Index along `v`, in `[0, degree()]`.
         * @return Const reference to that control point.
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
            std::vector<TPointT> reduced_along_v(nb_control_points());
            for (std::size_t i = 0; i < nb_control_points(); ++i) {
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

        /**
         * @brief Splits the surface along `u = AU` into the 2 surfaces of the same degree covering
         * `u` in `[0, AU]` and `[AU, 1]`.
         *
         * Exact, like the `BezierCurve::split()` it is built from: a tensor-product surface is a
         * Bezier curve in `u` whose control points are themselves curves in `v`, so subdividing every
         * `u`-fiber of the control grid subdivides the surface —
         * `first.value(u, v) == value(u * AU, v)` up to round-off, and the 2 halves share their
         * junction isoparametric curve exactly.
         *
         * @param AU Split parameter, strictly inside `(0, 1)` (see `BezierCurve::split()`).
         * @return The `[0, AU]` half first, the `[AU, 1]` half second.
         */
        std::pair<BezierSurface, BezierSurface> split_u(double AU) const {
            const std::size_t n = nb_control_points();
            ControlGrid low = make_grid(degree());
            ControlGrid high = make_grid(degree());
            for (std::size_t j = 0; j < n; ++j) {
                std::vector<TPointT> fiber(n);
                for (std::size_t i = 0; i < n; ++i) {
                    fiber[i] = m_control_points[i][j];
                }
                const auto [l, h] = BezierCurve<TPointT>(std::move(fiber)).split(AU);
                for (std::size_t i = 0; i < n; ++i) {
                    low[i][j] = l[i];
                    high[i][j] = h[i];
                }
            }
            return {BezierSurface(std::move(low)), BezierSurface(std::move(high))};
        }

        /**
         * @brief Splits the surface along `v = AV` into the 2 surfaces of the same degree covering
         * `v` in `[0, AV]` and `[AV, 1]`.
         * @param AV Split parameter, strictly inside `(0, 1)` (see `BezierCurve::split()`).
         * @return The `[0, AV]` half first, the `[AV, 1]` half second.
         * @see split_u() for why this is exact rather than a refit.
         */
        std::pair<BezierSurface, BezierSurface> split_v(double AV) const {
            const std::size_t n = nb_control_points();
            ControlGrid low = make_grid(degree());
            ControlGrid high = make_grid(degree());
            for (std::size_t i = 0; i < n; ++i) {
                const auto [l, h] = BezierCurve<TPointT>(m_control_points[i]).split(AV);
                low[i] = l.control_points();
                high[i] = h.control_points();
            }
            return {BezierSurface(std::move(low)), BezierSurface(std::move(high))};
        }

    private:
        /** @brief The `(degree()+1) x (degree()+1)` grid storing the control points. */
        ControlGrid m_control_points;

        /**
         * @brief One De Casteljau reduction of a row of control points.
         * @param APoints The control points to reduce (taken by value, reduced in place).
         * @param AT The parameter to reduce at.
         * @return The single point the reduction converges to.
         */
        static TPointT de_casteljau(std::vector<TPointT> APoints, double AT) {
            const std::size_t n = APoints.size() - 1;
            for (std::size_t k = 1; k <= n; ++k) {
                for (std::size_t i = 0; i <= n - k; ++i) {
                    APoints[i] = APoints[i] + (APoints[i + 1] - APoints[i]) * AT;
                }
            }
            return APoints[0];
        }
    };
} // namespace gecko
