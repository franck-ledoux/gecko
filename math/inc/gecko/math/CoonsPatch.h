#pragma once

#include <cstddef>

#include <gecko/math/BezierSurface.h>
#include <gecko/math/BezierVolume.h>
#include <gecko/math/CurveSurfaceTraits.h>
#include <gecko/math/Point3d.h>
#include <gecko/math/Vector3d.h>

namespace gecko {

    /**
     * @brief Builds the surface (via `CurveSurfaceTraits<TEdgeCurve>`, `BezierSurface` for
     * `BezierCurve` boundaries today) whose control net exactly reproduces the bilinearly-blended
     * Coons patch of 4 degree-`TEdgeCurve::Degree` boundary curves.
     *
     * This is the discrete (control-point-space) analog of the classic continuous Coons formula
     * `S(u,v) = (1-v)*Cu0(u) + v*Cu1(u) + (1-u)*Cv0(v) + u*Cv1(v) - bilinear(corners)`: applying the
     * same *ordinary* (unweighted) blend to each boundary curve's `i`-th/`j`-th *control point*
     * (instead of a continuous evaluation at `u=i/N`, `v=j/N`) yields the exact control net of the
     * surface a continuous Coons evaluation would trace — a standard CAGD result for polynomial
     * boundary curves, not an approximation. For `TEdgeCurve::Degree == 1` (straight edges) this
     * collapses exactly to bilinear corner interpolation. The return type is looked up generically
     * through `CurveSurfaceTraits`, so this function itself never names `BezierSurface`; the blend
     * below is still an unweighted control-point average, appropriate for `BezierCurve`/a uniform
     * B-spline boundary — a genuinely rational (NURBS) boundary would need a weighted variant of this
     * same construction, not just a new `CurveSurfaceTraits` specialization.
     *
     * @tparam TEdgeCurve Boundary curve representation (e.g. `BezierCurve<N, Point3d>`); must expose
     * `Degree`, an ordered, endpoint-inclusive `control_points()` (see `EdgeCurveConcept` in
     * `geom_itf`, not depended on here to keep `math` dependency-free of `geom_itf`), and a
     * `CurveSurfaceTraits<TEdgeCurve>::Surface` specialization.
     * @param AEdgeU0 Boundary curve along `u`, at `v=0` (from the `(u=0,v=0)` corner to `(u=1,v=0)`).
     * @param AEdgeU1 Boundary curve along `u`, at `v=1` (from `(u=0,v=1)` to `(u=1,v=1)`).
     * @param AEdgeV0 Boundary curve along `v`, at `u=0` (from `(u=0,v=0)` to `(u=0,v=1)`).
     * @param AEdgeV1 Boundary curve along `v`, at `u=1` (from `(u=1,v=0)` to `(u=1,v=1)`).
     * @return The surface reproducing the Coons blend of the 4 boundary curves.
     * @pre The 4 curves' shared corners must be the literal same point: `AEdgeU0`'s first/last
     * control points equal `AEdgeV0`'s/`AEdgeV1`'s first control points, and likewise for `AEdgeU1`.
     */
    template<typename TEdgeCurve>
    typename CurveSurfaceTraits<TEdgeCurve>::Surface coons_surface_from_edges(const TEdgeCurve &AEdgeU0,
                                                                              const TEdgeCurve &AEdgeU1,
                                                                              const TEdgeCurve &AEdgeV0,
                                                                              const TEdgeCurve &AEdgeV1) {
        constexpr std::size_t N = TEdgeCurve::Degree;
        constexpr std::size_t NumPts = N + 1;
        using Surface = typename CurveSurfaceTraits<TEdgeCurve>::Surface;

        const Point3d &origin = AEdgeU0.control_points()[0];
        typename Surface::ControlGrid grid;

        for (std::size_t i = 0; i < NumPts; ++i) {
            const double u = static_cast<double>(i) / static_cast<double>(N == 0 ? 1 : N);
            for (std::size_t j = 0; j < NumPts; ++j) {
                const double v = static_cast<double>(j) / static_cast<double>(N == 0 ? 1 : N);

                Vector3d acc(0.0, 0.0, 0.0);
                acc += Vector3d(origin, AEdgeU0.control_points()[i]) * (1.0 - v);
                acc += Vector3d(origin, AEdgeU1.control_points()[i]) * v;
                acc += Vector3d(origin, AEdgeV0.control_points()[j]) * (1.0 - u);
                acc += Vector3d(origin, AEdgeV1.control_points()[j]) * u;
                acc -= Vector3d(origin, AEdgeU0.control_points()[N]) * (u * (1.0 - v));
                acc -= Vector3d(origin, AEdgeU1.control_points()[0]) * ((1.0 - u) * v);
                acc -= Vector3d(origin, AEdgeU1.control_points()[N]) * (u * v);

                grid[i][j] = origin + acc;
            }
        }
        return Surface(grid);
    }

    /**
     * @brief Builds the volume (via `SurfaceVolumeTraits<TFaceSurface>`, `BezierVolume` for
     * `BezierSurface` faces today) whose control net exactly reproduces the 3D transfinite
     * interpolation ("TFI", sometimes loosely called a Coons volume) of its 6 bounding surfaces.
     *
     * Discrete analog of the classic Boolean-sum TFI construction `U + V + W - UV - VW - UW + UVW`
     * (`U`/`V`/`W` each ruled-blend a pair of opposite faces; `UV`/`VW`/`UW` are the double-counted
     * shared-edge contributions, subtracted; `UVW` is the triple-counted shared-corner contribution,
     * added back) — applied to each face's control point instead of a continuous evaluation, which
     * yields the exact control net (not an approximation) for boundary surfaces built the same
     * ordinary (unweighted) way as `coons_surface_from_edges()` — see that function's own note on
     * what a genuinely rational (NURBS) boundary would additionally require. Since each face is
     * itself built by coons_surface_from_edges(), its own boundary rows/columns already carry the
     * volume's 12 edges and 8 corners exactly, so only the 6 face grids are needed here —
     * `AFaceU0`/`AFaceU1` alone expose all 8 corners. For `TEdgeCurve::Degree == 1` this collapses
     * exactly to trilinear corner interpolation. The return type is looked up generically through
     * `SurfaceVolumeTraits`, so this function itself never names `BezierVolume`.
     *
     * @tparam TFaceSurface Face representation (e.g. `BezierSurface<N, Point3d>` as returned by
     * coons_surface_from_edges()); must expose `Degree`, a `(Degree+1) x (Degree+1)`
     * `control_points()` grid, and a `SurfaceVolumeTraits<TFaceSurface>::Volume` specialization.
     * @param AFaceU0 Bounding face at `u=0`, grid indexed `[v][w]`.
     * @param AFaceU1 Bounding face at `u=1`, grid indexed `[v][w]`.
     * @param AFaceV0 Bounding face at `v=0`, grid indexed `[u][w]`.
     * @param AFaceV1 Bounding face at `v=1`, grid indexed `[u][w]`.
     * @param AFaceW0 Bounding face at `w=0`, grid indexed `[u][v]`.
     * @param AFaceW1 Bounding face at `w=1`, grid indexed `[u][v]`.
     * @return The volume reproducing the 3D TFI blend of the 6 bounding faces.
     * @pre The 6 faces' shared boundary rows/columns (edges) and corners must be the literal same
     * points (guaranteed if all 6 were built by coons_surface_from_edges() from a shared set of 12
     * edge curves, each used by exactly 2 of the 6 faces).
     */
    template<typename TFaceSurface>
    typename SurfaceVolumeTraits<TFaceSurface>::Volume tfi_volume_from_faces(const TFaceSurface &AFaceU0,
                                                                             const TFaceSurface &AFaceU1,
                                                                             const TFaceSurface &AFaceV0,
                                                                             const TFaceSurface &AFaceV1,
                                                                             const TFaceSurface &AFaceW0,
                                                                             const TFaceSurface &AFaceW1) {
        constexpr std::size_t N = TFaceSurface::Degree;
        constexpr std::size_t NumPts = N + 1;
        using Volume = typename SurfaceVolumeTraits<TFaceSurface>::Volume;

        const Point3d &origin = AFaceU0.control_points()[0][0];
        typename Volume::ControlGrid grid;

        for (std::size_t i = 0; i < NumPts; ++i) {
            const double u = static_cast<double>(i) / static_cast<double>(N == 0 ? 1 : N);
            for (std::size_t j = 0; j < NumPts; ++j) {
                const double v = static_cast<double>(j) / static_cast<double>(N == 0 ? 1 : N);
                for (std::size_t k = 0; k < NumPts; ++k) {
                    const double w = static_cast<double>(k) / static_cast<double>(N == 0 ? 1 : N);

                    Vector3d acc(0.0, 0.0, 0.0);
                    const auto add = [&](const Point3d &p, double weight) { acc += Vector3d(origin, p) * weight; };

                    // Face (ruled) terms: U + V + W.
                    add(AFaceU0.control_points()[j][k], 1.0 - u);
                    add(AFaceU1.control_points()[j][k], u);
                    add(AFaceV0.control_points()[i][k], 1.0 - v);
                    add(AFaceV1.control_points()[i][k], v);
                    add(AFaceW0.control_points()[i][j], 1.0 - w);
                    add(AFaceW1.control_points()[i][j], w);

                    // Edge (double-counted) terms, subtracted: UV + VW + UW.
                    add(AFaceU0.control_points()[0][k], -(1.0 - u) * (1.0 - v));
                    add(AFaceU0.control_points()[N][k], -(1.0 - u) * v);
                    add(AFaceU1.control_points()[0][k], -u * (1.0 - v));
                    add(AFaceU1.control_points()[N][k], -u * v);

                    add(AFaceW0.control_points()[i][0], -(1.0 - v) * (1.0 - w));
                    add(AFaceW1.control_points()[i][0], -(1.0 - v) * w);
                    add(AFaceW0.control_points()[i][N], -v * (1.0 - w));
                    add(AFaceW1.control_points()[i][N], -v * w);

                    add(AFaceU0.control_points()[j][0], -(1.0 - u) * (1.0 - w));
                    add(AFaceU0.control_points()[j][N], -(1.0 - u) * w);
                    add(AFaceU1.control_points()[j][0], -u * (1.0 - w));
                    add(AFaceU1.control_points()[j][N], -u * w);

                    // Corner (triple-counted) term, added back: UVW.
                    add(AFaceU0.control_points()[0][0], (1.0 - u) * (1.0 - v) * (1.0 - w));
                    add(AFaceU0.control_points()[0][N], (1.0 - u) * (1.0 - v) * w);
                    add(AFaceU0.control_points()[N][0], (1.0 - u) * v * (1.0 - w));
                    add(AFaceU0.control_points()[N][N], (1.0 - u) * v * w);
                    add(AFaceU1.control_points()[0][0], u * (1.0 - v) * (1.0 - w));
                    add(AFaceU1.control_points()[0][N], u * (1.0 - v) * w);
                    add(AFaceU1.control_points()[N][0], u * v * (1.0 - w));
                    add(AFaceU1.control_points()[N][N], u * v * w);

                    grid[i][j][k] = origin + acc;
                }
            }
        }
        return Volume(grid);
    }

} // namespace gecko
