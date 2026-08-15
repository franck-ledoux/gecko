#pragma once
#include <algorithm>

#include <gecko/math/Point3d.h>
#include <gecko/math/Vector3d.h>

namespace gecko {

    /**
     * @brief Computes the closest point of a segment [@p a, @p b] to a query point.
     *
     * Standard clamped-projection formula: projects @p p onto the line through @p a and @p b,
     * then clamps the resulting parameter to [0, 1] so the result always lies on the segment.
     *
     * @param p Query point.
     * @param a First endpoint of the segment.
     * @param b Second endpoint of the segment.
     * @return The closest point of the segment to @p p.
     * @note If @p a and @p b coincide (a degenerate, zero-length segment), @p a is returned.
     */
    [[nodiscard]] inline Point3d
    closest_point_on_segment(const Point3d &p, const Point3d &a, const Point3d &b) noexcept {
        const Vector3d ab(a, b);
        const double denom = ab.norm_sq();
        if (denom <= 0.0) {
            return a; // Degenerate segment: a == b.
        }
        const Vector3d ap(a, p);
        const double t = std::clamp(dot(ap, ab) / denom, 0.0, 1.0);
        return a + t * ab;
    }

    /**
     * @brief Computes the closest point of a triangle (@p a, @p b, @p c) to a query point.
     *
     * Implements the barycentric/Voronoi-region method described in Christer Ericson's
     * "Real-Time Collision Detection" (section 5.1.5, "Closest Point on Triangle to Point"):
     * the point is classified into one of the triangle's 7 Voronoi regions (3 vertices, 3 edges,
     * 1 face) using signed areas derived from dot products, without computing the triangle's
     * plane or normal explicitly.
     *
     * @param p Query point.
     * @param a First corner of the triangle.
     * @param b Second corner of the triangle.
     * @param c Third corner of the triangle.
     * @return The closest point of the triangle to @p p.
     * @note If the triangle is degenerate (zero area, e.g. collinear corners), the closest point
     *       among its 3 edges (via closest_point_on_segment()) is returned instead.
     */
    [[nodiscard]] inline Point3d
    closest_point_on_triangle(const Point3d &p, const Point3d &a, const Point3d &b, const Point3d &c) noexcept {
        const Vector3d ab(a, b);
        const Vector3d ac(a, c);
        const Vector3d ap(a, p);
        const double d1 = dot(ab, ap);
        const double d2 = dot(ac, ap);
        if (d1 <= 0.0 && d2 <= 0.0) {
            return a; // Vertex region A.
        }

        const Vector3d bp(b, p);
        const double d3 = dot(ab, bp);
        const double d4 = dot(ac, bp);
        if (d3 >= 0.0 && d4 <= d3) {
            return b; // Vertex region B.
        }

        const double vc = d1 * d4 - d3 * d2;
        if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
            const double v = d1 / (d1 - d3);
            return a + v * ab; // Edge region AB.
        }

        const Vector3d cp(c, p);
        const double d5 = dot(ab, cp);
        const double d6 = dot(ac, cp);
        if (d6 >= 0.0 && d5 <= d6) {
            return c; // Vertex region C.
        }

        const double vb = d5 * d2 - d1 * d6;
        if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
            const double w = d2 / (d2 - d6);
            return a + w * ac; // Edge region AC.
        }

        const double va = d3 * d6 - d5 * d4;
        if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
            const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return b + w * Vector3d(b, c); // Edge region BC.
        }

        const double denom = va + vb + vc;
        if (denom <= 0.0) {
            // Degenerate (near-zero-area) triangle: fall back to the closest of its 3 edges.
            Point3d best = closest_point_on_segment(p, a, b);
            double best_dist = Vector3d(p, best).norm_sq();
            if (const Point3d on_bc = closest_point_on_segment(p, b, c); Vector3d(p, on_bc).norm_sq() < best_dist) {
                best = on_bc;
                best_dist = Vector3d(p, best).norm_sq();
            }
            if (const Point3d on_ca = closest_point_on_segment(p, c, a); Vector3d(p, on_ca).norm_sq() < best_dist) {
                best = on_ca;
            }
            return best;
        }

        // Face region: p projects inside the triangle.
        const double inv = 1.0 / denom;
        const double v = vb * inv;
        const double w = vc * inv;
        return a + v * ab + w * ac;
    }

} // namespace gecko
