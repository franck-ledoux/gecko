#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace gecko {
    /**
     * @class BezierCurve
     * @brief A Bezier curve of any degree, over control points of type @p TPointT.
     *
     * The degree is a *property of the object*, not of its type: one `BezierCurve<Point3d>` may be a
     * straight segment and the next a cubic, and either can be replaced by the other in place. That
     * is what lets a whole block structure change order while it is being edited — see
     * `Blocking::set_degree()` — which a degree baked into the type could not do without rebuilding
     * every cell into a different C++ type.
     *
     * @tparam TPointT The geometric point type, following gecko's point/vector affine-space split
     *                (e.g. gecko::Point3d): `TPointT - TPointT` yields a difference/vector type `V`,
     *                `V * double` yields `V`, and `TPointT + V` yields `TPointT`.
     * @note derivative()/value_derivative()/project()/projectPoint() are not written to that
     *       affine-space convention (they still assume `TPointT` itself supports `operator+`/
     *       `operator*(double)` directly, which gecko::Point3d does not) and are unused by the
     *       `block` module; only value()/operator()/split() are relied upon there.
     */
    template<typename TPointT>
    class BezierCurve {
    public:
        /**
         * @brief Default constructor: a degree-1 curve with both control points default-initialized.
         * @note Needed so this type can be used as a CGAL `Cell_attribute` info payload (default-
         *       constructible is a hard CGAL requirement); the curve is meaningless until its degree
         *       and control points are set.
         */
        BezierCurve() : m_control_points(2) {}

        /**
         * @brief Constructs a curve of a given degree, control points default-initialized.
         * @param ADegree The degree; the curve gets `ADegree + 1` control points.
         */
        explicit BezierCurve(std::size_t ADegree) : m_control_points(ADegree + 1) {}

        /**
         * @brief Constructs a curve from its control points; the degree follows from how many.
         * @param AControlPoints At least 2 control points.
         */
        explicit BezierCurve(std::vector<TPointT> AControlPoints) : m_control_points(std::move(AControlPoints)) {
            assert(m_control_points.size() >= 2 && "BezierCurve: a curve needs at least 2 control points");
        }

        /**
         * @brief Constructs a curve from its control points listed inline.
         * @param AControlPoints At least 2 control points.
         */
        BezierCurve(std::initializer_list<TPointT> AControlPoints) : m_control_points(AControlPoints) {
            assert(m_control_points.size() >= 2 && "BezierCurve: a curve needs at least 2 control points");
        }

        /**
         * @brief Constructs a curve from its control points passed as individual arguments.
         * @tparam TArgs Pack of @p TPointT, at least 2 of them.
         * @param AArgs The control points.
         */
        template<typename... TArgs>
            requires(sizeof...(TArgs) >= 2) && (std::same_as<std::decay_t<TArgs>, TPointT> && ...)
        explicit BezierCurve(TArgs &&...AArgs) : m_control_points{std::forward<TArgs>(AArgs)...} {}

        /** @brief Gets the curve's degree. @return The degree (control point count minus 1). */
        [[nodiscard]] std::size_t degree() const { return m_control_points.size() - 1; }

        /** @brief Gets how many control points define the curve. @return The count, `degree() + 1`. */
        [[nodiscard]] std::size_t nb_control_points() const { return m_control_points.size(); }

        /**
         * @brief Gets a read-only reference to all control points.
         * @return Const reference to the control points.
         */
        const std::vector<TPointT> &control_points() const { return m_control_points; }

        /**
         * @brief Accesses a control point by index.
         * @param AIndex Index of the control point, in `[0, degree()]`.
         * @return Reference to the requested control point.
         */
        TPointT &operator[](std::size_t AIndex) { return m_control_points[AIndex]; }

        /**
         * @brief Accesses a control point by index.
         * @param AIndex Index of the control point, in `[0, degree()]`.
         * @return Const reference to the requested control point.
         */
        const TPointT &operator[](std::size_t AIndex) const { return m_control_points[AIndex]; }

        /**
         * @brief Evaluates the point on the curve at parameter t using De Casteljau's algorithm.
         * De Casteljau's algorithm is numerically stable and avoids direct calculation
         * of Bernstein polynomials.
         *
         * @param AC curvilinear coordinate in the range [0.0, 1.0].
         * @return The evaluated point of type TPointT at parameter t.
         * @pre Parameter AC must satisfy 0.0 <= AC <= 1.0 (enforced via assert in Debug builds).
         */
        TPointT value(double AC) const {
            assert(AC >= 0.0 && AC <= 1.0 && "BezierCurve::value: Parameter AC must be in the range [0.0, 1.0]");
            std::vector<TPointT> points = m_control_points;
            const std::size_t n = degree();
            for (std::size_t k = 1; k <= n; ++k) {
                for (std::size_t i = 0; i <= n - k; ++i) {
                    points[i] = lerp(points[i], points[i + 1], AC);
                }
            }
            return points[0];
        }

        /**
         * @brief Evaluates the curve at parameter t (functor operator shortcut).
         * @param t Evaluation parameter, typically in the range [0.0, 1.0].
         * @return The evaluated point of type TPointT at parameter t.
         * @see value()
         */
        TPointT operator()(double t) const { return value(t); }

        /**
         * @brief Splits the curve at parameter @p AT into the 2 curves of the same degree covering
         * `[0, AT]` and `[AT, 1]`.
         *
         * This is De Casteljau *subdivision*, and it is exact rather than a fit: the 2 halves
         * reproduce the parent curve point for point, only reparameterized —
         * `first.value(t) == value(t * AT)` and `second.value(t) == value(AT + t * (1 - AT))`, up to
         * floating-point round-off. That exactness is what lets a block-structure cut keep the
         * geometry it cut through (see `Blocking::cut_sheet()`): re-deriving the halves from their
         * new boundaries instead would move the surface, since the restriction of a Coons patch to a
         * sub-rectangle is not the Coons patch of that sub-rectangle's own boundary curves.
         *
         * The 2 halves share their junction control point exactly (`first[degree()] == second[0]`),
         * which is the point at @p AT — so a caller placing a new node there can take it from either.
         *
         * The construction is the same De Casteljau pyramid `value()` runs, reading off its 2 outer
         * edges instead of only its apex: level `k`'s first point is the left half's `k`-th control
         * point, its last point the right half's.
         *
         * @param AT Split parameter, strictly inside `(0, 1)` — an endpoint would make one half
         *        degenerate, and callers are expected to reject that themselves with a real error
         *        rather than rely on this (enforced via assert in Debug builds only).
         * @return The `[0, AT]` half first, the `[AT, 1]` half second.
         */
        std::pair<BezierCurve, BezierCurve> split(double AT) const {
            assert(AT > 0.0 && AT < 1.0 && "BezierCurve::split: AT must be strictly inside (0.0, 1.0)");

            const std::size_t n = degree();
            std::vector<TPointT> points = m_control_points;
            std::vector<TPointT> left(n + 1);
            std::vector<TPointT> right(n + 1);
            left[0] = points[0];
            right[n] = points[n];

            for (std::size_t k = 1; k <= n; ++k) {
                for (std::size_t i = 0; i <= n - k; ++i) {
                    points[i] = lerp(points[i], points[i + 1], AT);
                }
                left[k] = points[0];
                right[n - k] = points[n - k];
            }
            return {BezierCurve(std::move(left)), BezierCurve(std::move(right))};
        }

        /**
         * @brief Computes the K-th derivative curve of the current Bezier curve.
         * @param AK The order of differentiation (e.g. 1 for first derivative, 2 for second).
         * @return A BezierCurve of degree `degree() - AK` representing the K-th derivative curve.
         * @pre AK must be less than or equal to degree().
         */
        [[nodiscard]] BezierCurve derivative(std::size_t AK = 1) const {
            assert(AK <= degree() && "BezierCurve::derivative: order higher than the curve degree");
            BezierCurve result = *this;
            for (std::size_t step = 0; step < AK; ++step) {
                const std::size_t n = result.degree();
                std::vector<TPointT> points(n);
                for (std::size_t i = 0; i < n; ++i) {
                    points[i] = (result.m_control_points[i + 1] - result.m_control_points[i]) * static_cast<double>(n);
                }
                result.m_control_points = std::move(points);
            }
            return result;
        }

        /**
         * @brief Evaluates the K-th derivative at parameter AT.
         * @param AT Evaluation parameter within [0.0, 1.0].
         * @param AK The order of differentiation (1 = velocity, 2 = acceleration, etc.).
         * @return The evaluated K-th derivative vector of type TPointT at parameter AT.
         */
        TPointT value_derivative(double AT, std::size_t AK = 1) const { return derivative(AK).value(AT); }

        /**
         * @brief Projects a given point P onto the Bezier curve.
         *
         * Finds the parameter t in [0.0, 1.0] that minimizes the distance between P and C(t).
         * Uses a combined approach: uniform sampling for initial guess, followed by
         * Newton-Raphson iterations to solve (C(t) - P) . C'(t) = 0.
         *
         * @param AP The query point P to project.
         * @param num_samples Number of initial sampling points to find a good starting point (default: 20).
         * @param max_iterations Maximum number of Newton-Raphson refinement steps (default: 10).
         * @param tolerance Convergence threshold on parameter t change (default: 1e-6).
         * @return Parameter t in range [0.0, 1.0] corresponding to the projected point C(t).
         */
        double project(const TPointT &AP,
                       std::size_t num_samples = 20,
                       std::size_t max_iterations = 10,
                       double tolerance = 1e-6) const {
            const BezierCurve d1 = derivative(1);
            const bool has_d2 = degree() >= 2;
            const BezierCurve d2 = has_d2 ? derivative(2) : BezierCurve();

            // 1. Initial Coarse Search: Sample t uniformly in [0, 1] to find the closest sample
            double best_t = 0.0;
            double min_sq_dist = std::numeric_limits<double>::max();

            for (std::size_t i = 0; i <= num_samples; ++i) {
                double t = static_cast<double>(i) / static_cast<double>(num_samples);
                TPointT c_t = value(t);
                TPointT diff = c_t - AP;
                double sq_dist = dot(diff, diff);

                if (sq_dist < min_sq_dist) {
                    min_sq_dist = sq_dist;
                    best_t = t;
                }
            }

            // 2. Refinement using Newton-Raphson on f(t) = (C(t) - P) . C'(t) = 0
            // f'(t) = ||C'(t)||^2 + (C(t) - P) . C''(t)
            double t = best_t;

            for (std::size_t iter = 0; iter < max_iterations; ++iter) {
                TPointT c = value(t);
                TPointT c_prime = d1.value(t);
                TPointT c_prime_prime = has_d2 ? d2.value(t) : TPointT{};
                TPointT diff = c - AP;

                double f = dot(diff, c_prime);
                double f_prime = dot(c_prime, c_prime) + dot(diff, c_prime_prime);

                if (std::abs(f_prime) < 1e-12) {
                    break;
                }

                double delta_t = f / f_prime;
                t -= delta_t;
                t = std::clamp(t, 0.0, 1.0);

                if (std::abs(delta_t) < tolerance) {
                    break; // Converged
                }
            }

            // 3. Final safety check: compare boundary endpoints t=0 and t=1 against the refined result
            double sq_dist_t = dot(value(t) - AP, value(t) - AP);
            double sq_dist_0 = dot(value(0.0) - AP, value(0.0) - AP);
            double sq_dist_1 = dot(value(1.0) - AP, value(1.0) - AP);

            if (sq_dist_0 < sq_dist_t && sq_dist_0 < sq_dist_1) return 0.0;
            if (sq_dist_1 < sq_dist_t) return 1.0;

            return t;
        }

        /**
         * @brief Projects a given point P and returns the projected TPointT on the curve.
         * @param point The query point P to project.
         * @param num_samples Number of initial sampling points to find a good starting point (default: 20).
         * @param max_iterations Maximum number of Newton-Raphson refinement steps (default: 10).
         * @param tolerance Convergence threshold on parameter t change (default: 1e-6).
         * @return The projected point C(t) on the curve, closest to @p point.
         * @see project()
         */
        TPointT projectPoint(const TPointT &point,
                             std::size_t num_samples = 20,
                             std::size_t max_iterations = 10,
                             double tolerance = 1e-6) const {
            double t_proj = project(point, num_samples, max_iterations, tolerance);
            return value(t_proj);
        }

    private:
        /** @brief The `degree() + 1` control points. */
        std::vector<TPointT> m_control_points;

        /**
         * @brief Performs linear interpolation between two points.
         * @param AP1 Start point (t = 0.0).
         * @param AP2 End point (t = 1.0).
         * @param AT Interpolation factor.
         * @return Interpolated point: AP1 + t * (AP2 - AP1), written via the point/vector affine
         *         operations (`-`/`+`/vector`*double`) rather than a direct weighted sum of points,
         *         since gecko::Point3d (the typical TPointT) has no `operator+(Point3d)` or
         *         `operator*(double)` of its own — only affine point-vector arithmetic.
         */
        static TPointT lerp(const TPointT &AP1, const TPointT &AP2, double AT) { return AP1 + (AP2 - AP1) * AT; }

        /**
         * @brief Computes the dot product between two points/vectors.
         * Uses dot() method if available (e.g., Eigen), otherwise falls back to coordinate multiplication.
         */
        static double dot(const TPointT &a, const TPointT &b) {
            if constexpr (requires { a.dot(b); }) {
                return a.dot(b);
            } else {
                return a.x * b.x + a.y * b.y;
            }
        }
    };
} // namespace gecko
