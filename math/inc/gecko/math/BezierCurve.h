#pragma once

#include <array>
#include <utility>
#include <algorithm>
#include <cassert>

namespace gecko {
    /**
     * @class BezierCurve
     * @brief Represents a Bezier curve of degree @param TN with control points of type @param TPointT.
     *
     * This class provides a Bezier curve representation where the degree is given as a template parameter
     * and it can be used for any type of points.
     *
     * @tparam TN The degree of the Bezier curve (number of control points is N + 1).
     * @tparam TPointT The geometric point type, following gecko's point/vector affine-space split
     *                (e.g. gecko::Point3d): `TPointT - TPointT` yields a difference/vector type `V`,
     *                `V * double` yields `V`, and `TPointT + V` yields `TPointT`.
     * @note derivative()/value_derivative()/project()/projectPoint() are not yet updated to this
     *       affine-space convention (they still assume `TPointT` itself supports `operator+`/
     *       `operator*(double)` directly, which gecko::Point3d does not) and are unused by the
     *       `block` module; only value()/operator() are relied upon there.
     */
    template<std::size_t TN, typename TPointT>
    class BezierCurve {
    public:
        /** @brief The degree of the Bezier curve. */
        static constexpr std::size_t Degree = TN;

        /** @brief The total number of control points defining the curve (N + 1). */
        static constexpr std::size_t NumControlPoints = TN + 1;

        /**
         * @brief Default constructor. Control points default-initialize to `TPointT{}`.
         * @note Needed so this type can be used as a CGAL `Cell_attribute` info payload (default-
         *       constructible is a hard CGAL requirement); the curve is meaningless until its
         *       control points are set via operator[]/control_points().
         */
        BezierCurve() = default;

        /**
         * @brief Constructor. Constructs a Bezier curve from an array of control points.
         * @param AControlPoints An array containing exactly N + 1 control points
         */
        explicit BezierCurve(const std::array<TPointT, NumControlPoints> &AControlPoints)
            : m_control_points(AControlPoints) {}

        /**
         * @brief Variadic constructor allowing direct inline initialization of control points.         *
         * @tparam Args Pack of point types.
         * @param args The N + 1 control points passed as individual arguments.
         * @note This constructor is enabled only if the number of arguments matches NumControlPoints.
         */
        template<typename... Args, typename = std::enable_if_t<sizeof...(Args) == NumControlPoints>>
        explicit BezierCurve(Args &&...args) : m_control_points{std::forward<Args>(args)...} {}

        /**
         * @brief Gets a read-only reference to all control points.
         * @return Const reference to the std::array of control points.
         */
        const std::array<TPointT, NumControlPoints> &control_points() const { return m_control_points; }

        /**
         * @brief Accesses a control point by index.
         * @param AIndex index of the control point (from 0 to N).
         * @return Reference to the requested control point.
         */
        TPointT &operator[](std::size_t AIndex) { return m_control_points[AIndex]; }

        /**
         * @brief Accesses a control point by index.
         * @param AIndex index of the control point (from 0 to N).
         * @return Const reference to the requested control point.
         */
        const TPointT &operator[](std::size_t AIndex) const { return m_control_points[AIndex]; }

        /**
         * @brief Evaluates the point on the curve at parameter t using De Casteljau's algorithm.
         * De Casteljau's algorithm is numerically stable and avoids direct calculation
         * of Bernstein polynomials.
         *
         * @param AC curvilinear coordinate in the range [0.0, 1.0].
         * @return The evaluated point of type PointT at parameter t.
         * @pre Parameter AC must satisfy 0.0 <= t <= 1.0 (enforced via assert in Debug builds).
         */
        TPointT value(double AC) const {
            assert(AC >= 0.0 && AC <= 1.0 && "BezierCurve::evaluate: Parameter AC must be in the range [0.0, 1.0]");
            // Copy the array of control points to work on it
            std::array<TPointT, NumControlPoints> points = m_control_points;

            // Successive reductions (De Casteljau)
            for (std::size_t k = 1; k <= TN; ++k) {
                for (std::size_t i = 0; i <= TN - k; ++i) {
                    points[i] = lerp(points[i], points[i + 1], AC);
                }
            }
            return points[0];
        }

        /**
         * @brief Evaluates the curve at parameter t (functor operator shortcut).
         * @param t Evaluation parameter, typically in the range [0.0, 1.0].
         * @return The evaluated point of type PointT at parameter t.
         * @see evaluate()
         */
        TPointT operator()(double t) const { return value(t); }

        /**
         * @brief Computes the K-th derivative curve of the current Bezier curve.
         *
         * Applies the derivative operation recursively K times at compile-time.
         *
         * @tparam TK The order of differentiation (e.g., 1 for first derivative, 2 for second).
         * @return A BezierCurve<N - K, PointT> representing the K-th derivative curve.
         * @pre K must be less than or equal to N (K <= N).
         */
        template<std::size_t TK>
        auto derivative() const {
            static_assert(TK <= TN, "Cannot compute a derivative of order higher than the curve degree.");

            if constexpr (TK == 0) {
                return *this;
            } else {
                return compute_derivative<TK>();
            }
        }

        /**
         * @brief Evaluates the K-th derivative at parameter AT.
         *
         * @tparam TK The order of differentiation (1 = velocity, 2 = acceleration, etc.).
         * @param AT Evaluation parameter within [0.0, 1.0].
         * @return The evaluated K-th derivative vector of type PointT at parameter t.
         * @pre Parameter AT must satisfy 0.0 <= AT <= 1.0.
         * @pre Order TK must satisfy TK <= TN.
         */
        template<std::size_t TK>
        TPointT value_derivative(double AT) const {
            return derivative<TK>().evaluate(AT);
        }

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
         * @pre Degree N must be at least 1.
         */
        double project(const TPointT &AP,
                       std::size_t num_samples = 20,
                       std::size_t max_iterations = 10,
                       double tolerance = 1e-6) const {
            static_assert(TN > 0, "Cannot project a point onto a 0-degree (single point) curve.");

            // Precompute first and second derivative curves
            const auto d1 = derivative<1>();

            // Handle degree 1 (line segment) vs higher degrees for second derivative
            auto eval_d2 = [&](double t) -> TPointT {
                if constexpr (TN >= 2) {
                    return derivative<2>().evaluate(t);
                } else {
                    return TPointT{}; // Second derivative is zero for a line
                }
            };

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
                TPointT c_prime_prime = eval_d2(t);
                TPointT diff = c - AP;

                // f(t) = (C(t) - P) . C'(t)
                double f = dot(diff, c_prime);

                // f'(t) = C'(t) . C'(t) + (C(t) - P) . C''(t)
                double f_prime = dot(c_prime, c_prime) + dot(diff, c_prime_prime);

                // Avoid division by zero
                if (std::abs(f_prime) < 1e-12) {
                    break;
                }

                double delta_t = f / f_prime;
                t -= delta_t;

                // Clamp t to [0, 1] bounds
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
         * @brief Projects a given point P and returns the projected PointT on the curve.
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
        /** @brief The array storing the N + 1 control points. */
        std::array<TPointT, NumControlPoints> m_control_points;

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
         * @brief Helper function to compute the K-th derivative curve recursively.
         */
        template<std::size_t TK>
        auto compute_derivative() const {
            if constexpr (TK == 1) {
                std::array<TPointT, TN> deriv_points;
                for (std::size_t i = 0; i < TN; ++i) {
                    deriv_points[i] = (m_control_points[i + 1] - m_control_points[i]) * static_cast<double>(TN);
                }

                return BezierCurve<TN - 1, TPointT>(deriv_points);
            } else {
                return derivative().template derivative<TK - 1>();
            }
        }

        /**
         * @brief Computes the dot product between two points/vectors.
         * Uses dot() method if available (e.g., Eigen), otherwise falls back to coordinate multiplication.
         */
        static double dot(const TPointT &a, const TPointT &b) {
            // En C++17, interpolation/produit scalaire standard pour x et y (Point2D simple)
            if constexpr (requires { a.dot(b); }) {
                return a.dot(b);
            } else {
                return a.x * b.x + a.y * b.y; // Fonctionne pour les structures simples avec .x et .y
            }
        }
    };
} // namespace gecko
