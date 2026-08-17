#pragma once

#include <gecko/math/BezierCurve.h>
#include <gecko/math/BezierSurface.h>
#include <gecko/math/BezierVolume.h>

namespace gecko {

    /**
     * @brief Maps a curve representation type to its paired tensor-product surface representation
     * (same degree, same point type).
     *
     * The extensibility seam letting `block::Blocking` and `CoonsPatch.h`'s construction functions
     * stay generic over the curve representation instead of hardcoding a specific family: today only
     * `BezierCurve` is specialized below, but a future `NurbsCurve`/`BSplineCurve` plugs in by adding
     * its own specialization here, touching no other file (`Blocking`, `coons_surface_from_edges()`
     * etc. all look the paired surface type up through this trait, never naming `BezierSurface`
     * directly). Deliberately left undefined for any curve type without a specialization below, so
     * using an unpaired curve type is a compile error, not a silently-wrong default.
     * @tparam TEdgeCurve The curve representation type to look up.
     */
    template<typename TEdgeCurve>
    struct CurveSurfaceTraits;

    /** @brief `BezierCurve<N, TPointT>`'s paired surface representation is `BezierSurface<N, TPointT>`. */
    template<std::size_t TN, typename TPointT>
    struct CurveSurfaceTraits<BezierCurve<TN, TPointT>> {
        /** @brief The paired surface representation. */
        using Surface = BezierSurface<TN, TPointT>;
    };

    /**
     * @brief Maps a surface representation type to its paired tensor-product volume representation
     * (same degree, same point type) — the volume-side analog of `CurveSurfaceTraits`, used by
     * `block::Blocking` and `tfi_volume_from_faces()` to stay generic over the face/block geometry
     * representation.
     * @tparam TFaceSurface The surface representation type to look up.
     */
    template<typename TFaceSurface>
    struct SurfaceVolumeTraits;

    /** @brief `BezierSurface<N, TPointT>`'s paired volume representation is `BezierVolume<N, TPointT>`. */
    template<std::size_t TN, typename TPointT>
    struct SurfaceVolumeTraits<BezierSurface<TN, TPointT>> {
        /** @brief The paired volume representation. */
        using Volume = BezierVolume<TN, TPointT>;
    };

} // namespace gecko
