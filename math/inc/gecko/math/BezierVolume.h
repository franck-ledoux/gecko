#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <utility>

#include <gecko/math/BezierCurve.h>

namespace gecko {
    /**
     * @class BezierVolume
     * @brief Represents a tensor-product Bezier volume of degree @p TN x @p TN x @p TN, with a
     * `(TN+1) x (TN+1) x (TN+1)` grid of control points of type @p TPointT.
     *
     * Generalizes BezierCurve's De Casteljau evaluation to 3 parametric dimensions: evaluating at
     * `(u,v,w)` reduces each control-point "pencil" along `w`, then the resulting grid along `v`,
     * then the resulting row along `u`.
     *
     * @tparam TN The degree of the volume in each parametric direction (control grid is
     * `(TN+1) x (TN+1) x (TN+1)`).
     * @tparam TPointT The geometric point type, following gecko's point/vector affine-space split
     * (e.g. gecko::Point3d): `TPointT - TPointT` yields a difference/vector type `V`, `V * double`
     * yields `V`, and `TPointT + V` yields `TPointT` (see BezierCurve.h).
     */
    template<std::size_t TN, typename TPointT>
    class BezierVolume {
    public:
        /** @brief The degree of the volume in each parametric direction. */
        static constexpr std::size_t Degree = TN;

        /** @brief The number of control points along each parametric direction (N + 1). */
        static constexpr std::size_t NumControlPoints = TN + 1;

        /** @brief Row-major control point grid: `m_control_points[i][j][k]`, `i`/`j`/`k` along `u`/`v`/`w`. */
        using ControlGrid =
            std::array<std::array<std::array<TPointT, NumControlPoints>, NumControlPoints>, NumControlPoints>;

        /** @brief Default constructor. Every control point default-initializes to `TPointT{}`. */
        BezierVolume() = default;

        /**
         * @brief Constructor. Constructs a Bezier volume from an explicit control grid.
         * @param AControlPoints A `(TN+1)^3` grid of control points, `[i][j][k]` along `u`/`v`/`w`.
         */
        explicit BezierVolume(const ControlGrid &AControlPoints) : m_control_points(AControlPoints) {}

        /**
         * @brief Gets a read-only reference to the control grid.
         * @return Const reference to the `(TN+1)^3` control point grid.
         */
        const ControlGrid &control_points() const { return m_control_points; }

        /**
         * @brief Accesses a control point by its 3D grid index.
         * @param AI Index along `u` (0 to N).
         * @param AJ Index along `v` (0 to N).
         * @param AK Index along `w` (0 to N).
         * @return Reference to the requested control point.
         * @note Named `control_point`, not `operator()`, to avoid an unresolvable overload
         *       ambiguity with the `(double,double,double)` evaluation `operator()` below for
         *       calls with integer-literal arguments (equally convertible to `std::size_t` or
         *       `double`).
         */
        TPointT &control_point(std::size_t AI, std::size_t AJ, std::size_t AK) { return m_control_points[AI][AJ][AK]; }

        /**
         * @brief Accesses a control point by its 3D grid index.
         * @param AI Index along `u` (0 to N).
         * @param AJ Index along `v` (0 to N).
         * @param AK Index along `w` (0 to N).
         * @return Const reference to the requested control point.
         */
        const TPointT &control_point(std::size_t AI, std::size_t AJ, std::size_t AK) const {
            return m_control_points[AI][AJ][AK];
        }

        /**
         * @brief Evaluates the point in the volume at parameters (u, v, w) using De Casteljau's
         * algorithm applied first along `w`, then `v`, then `u`.
         * @param AU Parametric coordinate along `u`, in the range [0.0, 1.0].
         * @param AV Parametric coordinate along `v`, in the range [0.0, 1.0].
         * @param AW Parametric coordinate along `w`, in the range [0.0, 1.0].
         * @return The evaluated point in the volume.
         * @pre AU, AV and AW must each satisfy 0.0 <= x <= 1.0 (enforced via assert in Debug builds).
         */
        TPointT value(double AU, double AV, double AW) const {
            assert(AU >= 0.0 && AU <= 1.0 && "BezierVolume::value: AU must be in the range [0.0, 1.0]");
            assert(AV >= 0.0 && AV <= 1.0 && "BezierVolume::value: AV must be in the range [0.0, 1.0]");
            assert(AW >= 0.0 && AW <= 1.0 && "BezierVolume::value: AW must be in the range [0.0, 1.0]");
            std::array<TPointT, NumControlPoints> reduced_along_u;
            for (std::size_t i = 0; i < NumControlPoints; ++i) {
                std::array<TPointT, NumControlPoints> reduced_along_v;
                for (std::size_t j = 0; j < NumControlPoints; ++j) {
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
         * @return The evaluated point in the volume.
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
            ControlGrid low{};
            ControlGrid high{};
            ControlGrid self = m_control_points;
            for (std::size_t a = 0; a < NumControlPoints; ++a) {
                for (std::size_t b = 0; b < NumControlPoints; ++b) {
                    std::array<TPointT, NumControlPoints> fiber{};
                    for (std::size_t i = 0; i < NumControlPoints; ++i) {
                        fiber[i] = AAt(self, a, b, i);
                    }
                    const auto [l, h] = BezierCurve<TN, TPointT>(fiber).split(AT);
                    for (std::size_t i = 0; i < NumControlPoints; ++i) {
                        AAt(low, a, b, i) = l[i];
                        AAt(high, a, b, i) = h[i];
                    }
                }
            }
            return {BezierVolume(low), BezierVolume(high)};
        }

        /** @brief The `(TN+1)^3` grid storing the control points. */
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
