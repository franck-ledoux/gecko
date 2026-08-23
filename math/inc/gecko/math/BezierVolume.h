#pragma once

#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

#include <gecko/math/BezierCurve.h>

namespace gecko {
    /**
     * @class BezierVolume
     * @brief A tensor-product Bezier volume of any degree, with a `(degree+1)^3` grid of control
     * points of type @p TPointT.
     *
     * The degree is a property of the object rather than of its type, for the reason given on
     * `BezierCurve`: a block structure has to be able to change order in place.
     *
     * @tparam TPointT The geometric point type, following gecko's point/vector affine-space split
     *                (e.g. gecko::Point3d).
     */
    template<typename TPointT>
    class BezierVolume {
    public:
        /** @brief Row-major control point grid: `grid[i][j][k]`, `i`/`j`/`k` along `u`/`v`/`w`. */
        using ControlGrid = std::vector<std::vector<std::vector<TPointT>>>;

        /**
         * @brief Builds an empty control grid for a given degree, every point default-initialized.
         * @param ADegree The degree in each parametric direction.
         * @return A `(ADegree+1)^3` grid.
         */
        static ControlGrid make_grid(std::size_t ADegree) {
            return ControlGrid(ADegree + 1,
                               std::vector<std::vector<TPointT>>(ADegree + 1, std::vector<TPointT>(ADegree + 1)));
        }

        /** @brief Default constructor: a degree-1 volume with every control point default-initialized. */
        BezierVolume() : m_control_points(make_grid(1)) {}

        /**
         * @brief Constructs a volume of a given degree, control points default-initialized.
         * @param ADegree The degree in each parametric direction.
         */
        explicit BezierVolume(std::size_t ADegree) : m_control_points(make_grid(ADegree)) {}

        /**
         * @brief Constructs a Bezier volume from an explicit control grid; the degree follows from
         * its size.
         * @param AControlPoints A cubic grid of control points, `[i][j][k]` along `u`/`v`/`w`.
         */
        explicit BezierVolume(ControlGrid AControlPoints) : m_control_points(std::move(AControlPoints)) {
            assert(m_control_points.size() >= 2 && "BezierVolume: the control grid must be at least 2x2x2");
        }

        /** @brief Gets the volume's degree in each parametric direction. @return The degree. */
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
         * @param AK Index along `w`, in `[0, degree()]`.
         * @return Reference to that control point.
         */
        TPointT &control_point(std::size_t AI, std::size_t AJ, std::size_t AK) { return m_control_points[AI][AJ][AK]; }

        /**
         * @brief Accesses one control point.
         * @param AI Index along `u`, in `[0, degree()]`.
         * @param AJ Index along `v`, in `[0, degree()]`.
         * @param AK Index along `w`, in `[0, degree()]`.
         * @return Const reference to that control point.
         */
        const TPointT &control_point(std::size_t AI, std::size_t AJ, std::size_t AK) const {
            return m_control_points[AI][AJ][AK];
        }

        /**
         * @brief Evaluates the point on the volume at parameters (u, v, w) using De Casteljau's
         * algorithm applied along `w`, then `v`, then `u`.
         * @param AU Parametric coordinate along `u`, in the range [0.0, 1.0].
         * @param AV Parametric coordinate along `v`, in the range [0.0, 1.0].
         * @param AW Parametric coordinate along `w`, in the range [0.0, 1.0].
         * @return The evaluated point inside the volume.
         * @pre AU, AV and AW must each satisfy 0.0 <= x <= 1.0 (enforced via assert in Debug builds).
         */
        TPointT value(double AU, double AV, double AW) const {
            assert(AU >= 0.0 && AU <= 1.0 && "BezierVolume::value: AU must be in the range [0.0, 1.0]");
            assert(AV >= 0.0 && AV <= 1.0 && "BezierVolume::value: AV must be in the range [0.0, 1.0]");
            assert(AW >= 0.0 && AW <= 1.0 && "BezierVolume::value: AW must be in the range [0.0, 1.0]");
            const std::size_t n = nb_control_points();
            std::vector<TPointT> reduced_along_u(n);
            for (std::size_t i = 0; i < n; ++i) {
                std::vector<TPointT> reduced_along_v(n);
                for (std::size_t j = 0; j < n; ++j) {
                    reduced_along_v[j] = de_casteljau(m_control_points[i][j], AW);
                }
                reduced_along_u[i] = de_casteljau(reduced_along_v, AV);
            }
            return de_casteljau(reduced_along_u, AU);
        }

        /**
         * @brief Evaluates the volume at parameters (u, v, w) (functor operator shortcut).
         * @param AU Parametric coordinate along `u`.
         * @param AV Parametric coordinate along `v`.
         * @param AW Parametric coordinate along `w`.
         * @return The evaluated point inside the volume.
         * @see value()
         */
        TPointT operator()(double AU, double AV, double AW) const { return value(AU, AV, AW); }

        /**
         * @brief Splits the volume along `u = AU` into the 2 volumes of the same degree covering
         * `u` in `[0, AU]` and `[AU, 1]`.
         *
         * Exact, like the `BezierCurve::split()` it is built from: a tensor-product volume is a
         * Bezier curve in `u` whose control points are themselves surfaces in `(v, w)`, so
         * subdividing every `u`-fiber of the control grid subdivides the volume —
         * `first.value(u, v, w) == value(u * AU, v, w)` up to round-off, and the 2 halves share their
         * junction isoparametric surface exactly. That is what lets `Blocking::cut_sheet()` split a
         * curved block without moving the geometry it cut through.
         *
         * @param AU Split parameter, strictly inside `(0, 1)` (see `BezierCurve::split()`).
         * @return The `[0, AU]` half first, the `[AU, 1]` half second.
         */
        std::pair<BezierVolume, BezierVolume> split_u(double AU) const {
            return split_along(AU,
                               [](ControlGrid &AGrid, std::size_t AA, std::size_t AB, std::size_t AIdx) -> TPointT & {
                                   return AGrid[AIdx][AA][AB];
                               });
        }

        /**
         * @brief Splits the volume along `v = AV` into the 2 volumes of the same degree covering
         * `v` in `[0, AV]` and `[AV, 1]`.
         * @param AV Split parameter, strictly inside `(0, 1)` (see `BezierCurve::split()`).
         * @return The `[0, AV]` half first, the `[AV, 1]` half second.
         * @see split_u() for why this is exact rather than a refit.
         */
        std::pair<BezierVolume, BezierVolume> split_v(double AV) const {
            return split_along(AV,
                               [](ControlGrid &AGrid, std::size_t AA, std::size_t AB, std::size_t AIdx) -> TPointT & {
                                   return AGrid[AA][AIdx][AB];
                               });
        }

        /**
         * @brief Splits the volume along `w = AW` into the 2 volumes of the same degree covering
         * `w` in `[0, AW]` and `[AW, 1]`.
         * @param AW Split parameter, strictly inside `(0, 1)` (see `BezierCurve::split()`).
         * @return The `[0, AW]` half first, the `[AW, 1]` half second.
         * @see split_u() for why this is exact rather than a refit.
         */
        std::pair<BezierVolume, BezierVolume> split_w(double AW) const {
            return split_along(AW,
                               [](ControlGrid &AGrid, std::size_t AA, std::size_t AB, std::size_t AIdx) -> TPointT & {
                                   return AGrid[AA][AB][AIdx];
                               });
        }

    private:
        /**
         * @brief Shared body of `split_u()`/`split_v()`/`split_w()`: subdivides every control-grid
         * fiber running along one axis, the 3 axes differing only in how a fiber index maps back
         * into the grid.
         * @tparam TAt Callable `(ControlGrid&, std::size_t, std::size_t, std::size_t) -> TPointT&`
         *         addressing the grid entry at fiber position 3 on the fiber named by 1 and 2.
         * @param AT Split parameter, strictly inside `(0, 1)`.
         * @param AAt The grid-addressing callable.
         * @return The `[0, AT]` half first, the `[AT, 1]` half second.
         */
        template<typename TAt>
        std::pair<BezierVolume, BezierVolume> split_along(double AT, TAt AAt) const {
            const std::size_t n = nb_control_points();
            ControlGrid low = make_grid(degree());
            ControlGrid high = make_grid(degree());
            ControlGrid self = m_control_points;
            for (std::size_t a = 0; a < n; ++a) {
                for (std::size_t b = 0; b < n; ++b) {
                    std::vector<TPointT> fiber(n);
                    for (std::size_t i = 0; i < n; ++i) {
                        fiber[i] = AAt(self, a, b, i);
                    }
                    const auto [l, h] = BezierCurve<TPointT>(std::move(fiber)).split(AT);
                    for (std::size_t i = 0; i < n; ++i) {
                        AAt(low, a, b, i) = l[i];
                        AAt(high, a, b, i) = h[i];
                    }
                }
            }
            return {BezierVolume(std::move(low)), BezierVolume(std::move(high))};
        }

        /** @brief The `(degree()+1)^3` grid storing the control points. */
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
