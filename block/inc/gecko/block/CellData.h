#pragma once

#include <array>
#include <tuple>
#include <utility>
#include <vector>

#include <CGAL/Cell_attribute.h>
#include <CGAL/Combinatorial_map.h>
#include <gecko/math/Point3d.h>
#include <gecko/utils/Groups.h>

namespace gecko {

    /**
     * @struct CellInfo
     * @brief Classification data shared by every cell (node/edge/face/block) of a Blocking
     * structure: the set of geometric entities the cell may be classified on.
     *
     * Every entry's dimension is `>=` the cell's own topological dimension (a node may be
     * classified on a vertex, curve, surface or volume; an edge only on a curve, surface or
     * volume; never a lower-dimension entity, which would over-constrain — e.g. an edge is never
     * pinned to a single point). Empty = unclassified; size 1 is the common case; size > 1 when
     * several entities legitimately compete (e.g. a corner shared by several curves) — see
     * `Blocking::classify()`. Self-contained: no dependency on a geometric model's own `GroupId`/
     * `groups()`, so classification works the same whether or not the input model defines any
     * `$PhysicalNames` groups.
     */
    struct CellInfo {
        /** @brief The geometric entities this cell may be classified on, as (dimension, entity_tag) pairs. */
        std::vector<std::pair<GroupDim, Int>> geom_targets;
    };

    /** @brief Classification + position data for a Blocking node (0-cell). */
    struct NodeInfo : CellInfo {
        /** @brief Spatial location of the node. */
        Point3d point;
        /**
         * @brief A number identifying this node for as long as it exists, handed out by
         * `Blocking::create_node()` and never reused.
         *
         * There to give the class a *stable* total order on nodes, which the canonical local frames
         * of faces and blocks are built from (see `Blocking::face_local_nodes()`). Attribute handles
         * cannot serve: CGAL discards and rebuilds an attribute whenever the orbit behind it is
         * disturbed — which deleting a neighbouring block does — and a frame ordered by handles then
         * silently rotates underneath the surface stored in it. An id copied by `SplitFunctor`
         * survives exactly those rebuilds.
         */
        Int id = -1;
    };

    /**
     * @brief Classification + geometry data for a Blocking edge (1-cell).
     * @tparam TEdgeCurve Curve representation (see `EdgeCurveConcept` in `geom_itf`), e.g.
     * `BezierCurve<N, Point3d>`.
     */
    template<typename TEdgeCurve>
    struct EdgeInfo : CellInfo {
        /** @brief The edge's own curve, `N=1` (straight) by default until built/fitted. */
        TEdgeCurve curve;
    };

    /**
     * @brief Classification + geometry data for a Blocking face (2-cell): a genuine 3-cell's
     * bounding face in a 3D (hex) block, or the top-cell itself for a standalone 2D (quad) block.
     * @tparam TFaceSurface Surface representation, e.g. `BezierSurface<N, Point3d>`.
     */
    template<typename TFaceSurface>
    struct FaceInfo : CellInfo {
        /** @brief The face's own surface, built via Coons construction from its 4 boundary edges. */
        TFaceSurface surface;
        /**
         * @brief The `NodeInfo::id` of each of the 4 corners `surface` is parameterized against, at
         * `(u,v)` = (0,0), (1,0), (1,1), (0,1); all -1 until the face has a surface.
         *
         * Recorded rather than re-derived, and that is the whole point. A local frame worked out
         * from the face on demand — from its dart, from an ordering of its corners — rests on
         * something CGAL is free to change underneath it, and when it changes the surface stored in
         * that frame quietly reads back rotated or mirrored. Writing the frame down alongside the
         * surface it belongs to removes the question: whoever reads it later matches these ids
         * against the face's corners and gets the frame the surface was actually written in.
         */
        std::array<Int, 4> frame{-1, -1, -1, -1};
    };

    /**
     * @brief Classification + geometry data for a Blocking block (3-cell). Only meaningful for
     * genuine hex blocks — a standalone 2D (quad) block has no 3-cell at all.
     * @tparam TBlockVolume Volume representation, e.g. `BezierVolume<N, Point3d>`.
     */
    template<typename TBlockVolume>
    struct BlockInfo : CellInfo {
        /** @brief The block's own volume, built via TFI construction from its 6 bounding faces. */
        TBlockVolume volume;
        /** @brief The `NodeInfo::id` of each of the 8 corners `volume` is parameterized against, in
         * `HEX_CORNER_UVW` order; all -1 until the block has a volume.
         * @see FaceInfo::frame for why the frame is recorded and not re-derived. */
        std::array<Int, 8> frame{-1, -1, -1, -1, -1, -1, -1, -1};
    };

    /**
     * @struct MergeFunctor
     * @brief No-op `Cell_attribute` on-merge functor for edge/face/block attributes: keeps the
     * first attribute's data unchanged, discards the second.
     *
     * Safe because `Blocking::classify()` only ever runs *after* `build_connectivity()`'s auto-sew
     * (see `Blocking.h`), so no attribute ever being merged here already carries real classification
     * or geometry data — both sides are still in their freshly-created, default/straight state.
     */
    struct MergeFunctor {
        /**
         * @brief Merges 2 cell attributes by keeping the first one's data, in place, unchanged.
         * @tparam TCellAttribute A `CGAL::Cell_attribute` instantiation.
         */
        template<class TCellAttribute>
        void operator()(TCellAttribute & /*ACA1*/, TCellAttribute & /*ACA2*/) const {}
    };

    /**
     * @struct SplitFunctor
     * @brief `Cell_attribute` on-split functor for edge/face/block attributes: duplicates the
     * first attribute's data onto the second.
     *
     * Not exercised by any operation this module currently implements (no sheet-cut/pillow-style
     * topology edits), but CGAL's `Cell_attribute` requires a valid functor type regardless.
     */
    /**
     * @struct NodeSplitFunctor
     * @brief `Cell_attribute` on-split functor for node attributes: copies the original's data,
     * id included.
     *
     * CGAL splits a node attribute whenever the vertex orbit behind it comes apart — which deleting
     * a neighbouring block does — and the copy is still, geometrically and to every cell that has
     * one as a corner, the same corner. Keeping the id is what lets a cell find its own recorded
     * frame again afterwards (see `FaceInfo::frame`): the ids in a frame are matched against one
     * cell's own corners, where 2 copies of a single original node never appear together, so they
     * need only be unique within a cell and not across the map.
     */
    struct NodeSplitFunctor {
        /**
         * @brief Splits a node attribute, copying its data onto the new one but leaving that one's
         * id unset.
         * @tparam TCellAttribute A `CGAL::Cell_attribute` instantiation over `NodeInfo`.
         * @param ACA1 The original attribute, whose data is copied.
         * @param ACA2 The newly-created attribute.
         */
        template<class TCellAttribute>
        void operator()(TCellAttribute &ACA1, TCellAttribute &ACA2) const {
            ACA2.info() = ACA1.info();
        }
    };

    struct SplitFunctor {
        /**
         * @brief Splits a cell attribute by copying its data onto the newly-created second one.
         * @tparam TCellAttribute A `CGAL::Cell_attribute` instantiation.
         * @param ACA1 The original attribute, whose data is copied.
         * @param ACA2 The newly-created attribute, overwritten with a copy of @p ACA1's data.
         */
        template<class TCellAttribute>
        void operator()(TCellAttribute &ACA1, TCellAttribute &ACA2) const {
            ACA2.info() = ACA1.info();
        }
    };

    /**
     * @class CellData
     * @brief CGAL `Combinatorial_map` `Items` traits: describes the 4 kinds of cell attribute
     * (node/edge/face/block) a Blocking's underlying map stores.
     *
     * The map is always dimension 3 (see `Blocking.h`'s "one map, always 3D" design): a "2D" quad
     * block is simply a standalone face with no incident block attribute, not a different map
     * dimension, so `Attributes` unconditionally has all 4 entries.
     *
     * @tparam TEdgeCurve Edge curve representation, e.g. `BezierCurve<N, Point3d>`.
     * @tparam TFaceSurface Face surface representation, e.g. `BezierSurface<N, Point3d>`.
     * @tparam TBlockVolume Block volume representation, e.g. `BezierVolume<N, Point3d>`.
     */
    template<typename TEdgeCurve, typename TFaceSurface, typename TBlockVolume>
    struct CellData {
        /** @brief CGAL's expected nested attribute-declaration template. */
        template<class TCMap>
        struct Dart_wrapper {
            /** @brief Node (0-cell) attribute type. */
            using Node_attr = CGAL::Cell_attribute<TCMap, NodeInfo, CGAL::Tag_true, MergeFunctor, NodeSplitFunctor>;
            /** @brief Edge (1-cell) attribute type. */
            using Edge_attr =
                CGAL::Cell_attribute<TCMap, EdgeInfo<TEdgeCurve>, CGAL::Tag_true, MergeFunctor, SplitFunctor>;
            /** @brief Face (2-cell) attribute type. */
            using Face_attr =
                CGAL::Cell_attribute<TCMap, FaceInfo<TFaceSurface>, CGAL::Tag_true, MergeFunctor, SplitFunctor>;
            /** @brief Block (3-cell) attribute type. */
            using Block_attr =
                CGAL::Cell_attribute<TCMap, BlockInfo<TBlockVolume>, CGAL::Tag_true, MergeFunctor, SplitFunctor>;
            /** @brief Ordered tuple of attributes, one per dimension 0 to 3. */
            using Attributes = std::tuple<Node_attr, Edge_attr, Face_attr, Block_attr>;
        };
    };

    /**
     * @brief The (always dimension-3) combinatorial map type underlying a `Blocking<TGeomModel,
     * TEdgeCurve>` instance.
     * @tparam TEdgeCurve Edge curve representation, e.g. `BezierCurve<N, Point3d>`.
     * @tparam TFaceSurface Face surface representation, e.g. `BezierSurface<N, Point3d>`.
     * @tparam TBlockVolume Block volume representation, e.g. `BezierVolume<N, Point3d>`.
     */
    template<typename TEdgeCurve, typename TFaceSurface, typename TBlockVolume>
    using CMap = CGAL::Combinatorial_map<3, CellData<TEdgeCurve, TFaceSurface, TBlockVolume>>;

} // namespace gecko
