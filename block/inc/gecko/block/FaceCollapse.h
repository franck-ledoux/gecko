#pragma once

#include <utility>
#include <vector>

#include <gecko/block/Blocking.h>
#include <gecko/geom_itf/Concepts.h>
#include <gecko/math/Point3d.h>
#include <gecko/math/Vector3d.h>
#include <gecko/utils/Groups.h>

namespace gecko {

    /**
     * @class FaceCollapse
     * @brief 2D block-editing operation: collapses a quad face along one of its diagonals.
     *
     * The 2 diagonally-opposite selected corners merge into a single node; the 2 pairs of edges
     * that become duplicates as a result (each pair sharing the merged node and one of the quad's
     * other 2 corners) each merge into one edge too. An operation class, not a `Blocking` member
     * (per the data/operations split — see the design plan): it operates on `Blocking` purely
     * through its existing public API (`cmap()`, `geom_model()`, `rebuild_face_geometry()`,
     * `straight_curve()`, the `Face`/`Edge`/`Node` handle aliases).
     *
     * Classification: merging the 2 corners (and each merging edge pair) is only allowed when
     * their `geom_targets` are identical, or one side is unclassified (empty) — any other conflict
     * blocks the *whole* collapse via `can_apply()`, exactly like `ChordRemoval`'s coherence rule.
     *
     * Geometry: the merged node's position is the midpoint of the 2 collapsed corners, projected
     * onto the merged classification target if non-empty. Each surviving merged edge is rebuilt as
     * a straight line between its new endpoints (the 2 original edges' curves aren't blended — a
     * collapse fundamentally changes the quad's shape, so a fresh straight edge is simplest, matching
     * how newly-created edges elsewhere in this module always start straight too); any face whose
     * boundary includes such a merged edge has its surface rebuilt via `rebuild_face_geometry()`.
     *
     * V1 scope limit (documented, not silently skipped): `apply()` only refreshes the cells
     * *directly* touched by the collapse (the merged node, the 2 merged edges, their incident
     * faces) — it does not walk the whole structure to refresh other, more distant edges/faces
     * that merely share the moved node's old position. Call `Blocking::classify()` afterward for
     * full geometric consistency, exactly as it already is the designated pass for that.
     *
     * @tparam TGeomModel Geometric model type, must satisfy `GeomModelConcept`.
     * @tparam TEdgeCurve Edge curve representation, must satisfy `EdgeCurveConcept`.
     */
    template<GeomModelConcept TGeomModel, EdgeCurveConcept TEdgeCurve = BezierCurve<1, Point3d>>
    class FaceCollapse {
    public:
        /** @brief The `Blocking` specialization this operation acts on. */
        using BlockingT = Blocking<TGeomModel, TEdgeCurve>;
        /** @brief Handle to a face, as used by `BlockingT`. */
        using Face = typename BlockingT::Face;
        /** @brief Handle to an edge, as used by `BlockingT`. */
        using Edge = typename BlockingT::Edge;
        /** @brief Handle to a node, as used by `BlockingT`. */
        using Node = typename BlockingT::Node;
        /** @brief A dart descriptor, as used by `BlockingT`. */
        using Dart = typename BlockingT::Dart;

        /**
         * @brief Constructor.
         * @param ABlocking The blocking this operation will act on; must outlive this object.
         */
        explicit FaceCollapse(BlockingT &ABlocking) : m_blocking(ABlocking) {}

        /**
         * @brief Checks whether `AFace` can be collapsed along the diagonal `(ANode0, ANode2)`.
         *
         * Requires: the blocking is purely 2D; `ANode0` is one of `AFace`'s 4 corners and
         * `ANode2` is its diagonally-opposite corner (not an adjacent one); `ANode0`/`ANode2`'s
         * `geom_targets` are coherent (identical, or one empty); and each of the 2 resulting
         * merged-edge-pairs' `geom_targets` are coherent too.
         *
         * @param AFace The face to check.
         * @param ANode0 One corner of the diagonal.
         * @param ANode2 The other corner of the diagonal (must be diagonally opposite `ANode0`).
         * @return true if the collapse can be applied.
         */
        [[nodiscard]] bool can_apply(Face AFace, Node ANode0, Node ANode2) const {
            if (!m_blocking.is_purely_2d()) return false;
            auto &cmap = m_blocking.cmap();

            const Dart d0 = find_corner_dart(AFace, ANode0);
            if (d0 == cmap.null_descriptor) return false;
            const Dart d2 = cmap.template beta<1, 1>(d0);
            if (cmap.template attribute<0>(d2) != ANode2) return false;
            if (!coherent(ANode0->info().geom_targets, ANode2->info().geom_targets)) return false;

            const Dart d1 = cmap.template beta<1>(d0);
            const Dart d3 = cmap.template beta<0>(d0);
            const Edge e01 = cmap.template attribute<1>(d0);
            const Edge e12 = cmap.template attribute<1>(d1);
            const Edge e23 = cmap.template attribute<1>(d2);
            const Edge e30 = cmap.template attribute<1>(d3);
            if (!coherent(e01->info().geom_targets, e12->info().geom_targets)) return false;
            if (!coherent(e23->info().geom_targets, e30->info().geom_targets)) return false;
            return true;
        }

        /**
         * @brief Collapses `AFace` along the diagonal `(ANode0, ANode2)`.
         * @param AFace The face to collapse.
         * @param ANode0 One corner of the diagonal.
         * @param ANode2 The other corner of the diagonal (must be diagonally opposite `ANode0`).
         * @pre `can_apply(AFace, ANode0, ANode2)`
         */
        void apply(Face AFace, Node ANode0, Node ANode2) {
            assert(can_apply(AFace, ANode0, ANode2) && "FaceCollapse::apply: precondition violated");
            auto &cmap = m_blocking.cmap();

            const Dart d0 = find_corner_dart(AFace, ANode0);
            const Dart d1 = cmap.template beta<1>(d0);
            const Dart d2 = cmap.template beta<1, 1>(d0);
            const Dart d3 = cmap.template beta<0>(d0);

            // --- Precompute the coherent merge results before any topology mutation ---
            const auto merged_node_targets = coherent_merge(ANode0->info().geom_targets, ANode2->info().geom_targets);
            Point3d merged_pos = ANode0->info().point + Vector3d(ANode0->info().point, ANode2->info().point) * 0.5;
            if (!merged_node_targets.empty()) {
                merged_pos = project_onto_target(m_blocking.geom_model(),
                                                 merged_node_targets.front().first,
                                                 merged_node_targets.front().second,
                                                 merged_pos);
            }
            const Edge e01 = cmap.template attribute<1>(d0);
            const Edge e12 = cmap.template attribute<1>(d1);
            const Edge e23 = cmap.template attribute<1>(d2);
            const Edge e30 = cmap.template attribute<1>(d3);
            const auto merged_edge1_targets = coherent_merge(e01->info().geom_targets, e12->info().geom_targets);
            const auto merged_edge2_targets = coherent_merge(e23->info().geom_targets, e30->info().geom_targets);

            // --- Topology surgery: split the quad along the diagonal, then contract it ---
            const Dart diag_dart = cmap.insert_cell_1_in_cell_2(d0, d2);
            cmap.template contract_cell<1>(diag_dart);
            // d0/d1/d2/d3 remain valid dart handles (contract_cell<1> only erases the diagonal's
            // own 2 darts): {d0,d1} now bound one degenerate bigon face, {d2,d3} the other.

            // contract_cell<1> only propagates the vertex merge within ONE of the 2 triangles the
            // diagonal split the quad into (verified against CGAL's own Contract_cell_functor<1>:
            // its Group_attribute_functor_run<0,1> call only reconciles the *local* neighbor pair
            // adjacent to the contracted edge's own darts on the side it happens to process first —
            // it does not walk the wider vertex orbit to also fix up the *other* triangle's own
            // dart at the same nominal corner). d2 is exactly that other triangle's own dart at the
            // merged corner, so if it didn't already end up sharing d0's (whichever survived)
            // attribute, force it to explicitly — verified empirically with a standalone debug
            // program before trusting it, after this exact gap silently left 2 distinct node
            // attributes at the same position instead of one merged node.
            if (cmap.template attribute<0>(d0) != cmap.template attribute<0>(d2)) {
                cmap.template set_attribute<0>(d2, cmap.template attribute<0>(d0));
            }

            // Overwrite the merged node's info() right away, while d0 is still guaranteed valid —
            // collapse_bigon() below erases d0/d1/d2/d3 (and, in the fully-standalone case, may
            // garbage-collect the merged node attribute itself once nothing references it any
            // more), so there'd be no safe dart left to re-look it up afterward. Writing now is
            // harmless even if the node doesn't survive the bigon collapses: it simply gets freed
            // right after, with no dangling access ever happening in either case.
            cmap.template attribute<0>(d0)->info().point = merged_pos;
            cmap.template attribute<0>(d0)->info().geom_targets = merged_node_targets;

            std::vector<Face> faces_to_rebuild;
            const Edge survivor1 = collapse_bigon(d0, faces_to_rebuild);
            const Edge survivor2 = collapse_bigon(d2, faces_to_rebuild);

            // --- Overwrite surviving edges with the precomputed, coherent data ---
            if (survivor1 != cmap.null_descriptor) {
                survivor1->info().geom_targets = merged_edge1_targets;
                survivor1->info().curve = rebuilt_curve(survivor1, cmap);
            }
            if (survivor2 != cmap.null_descriptor) {
                survivor2->info().geom_targets = merged_edge2_targets;
                survivor2->info().curve = rebuilt_curve(survivor2, cmap);
            }
            for (const Face &f : faces_to_rebuild) {
                m_blocking.rebuild_face_geometry(f);
            }
        }

    private:
        /** @brief The `Blocking::Map` type, for brevity in helper signatures below. */
        using Map = typename BlockingT::Map;

        /**
         * @brief Finds the dart of `AFace`'s own boundary whose start corner is `ANode`.
         * @param AFace The face to walk.
         * @param ANode The corner to find.
         * @return The matching dart, or `cmap().null_descriptor` if `ANode` isn't one of
         *         `AFace`'s corners.
         */
        Dart find_corner_dart(Face AFace, Node ANode) const {
            auto &cmap = m_blocking.cmap();
            Dart walk = AFace->dart();
            for (int c = 0; c < 4; ++c) {
                if (cmap.template attribute<0>(walk) == ANode) return walk;
                walk = cmap.template beta<1>(walk);
            }
            return cmap.null_descriptor;
        }

        /**
         * @brief Checks whether 2 classification target sets are coherent: identical, or one
         * empty (adopting the other unconditionally).
         * @param A First target set.
         * @param B Second target set.
         * @return true if `A` and `B` don't conflict.
         */
        static bool coherent(const std::vector<std::pair<GroupDim, Int>> &A,
                             const std::vector<std::pair<GroupDim, Int>> &B) {
            return A.empty() || B.empty() || A == B;
        }

        /**
         * @brief Merges 2 already-`coherent()` classification target sets: whichever is
         * non-empty, or either (they're equal) if both are.
         * @param A First target set.
         * @param B Second target set.
         * @return The merged target set.
         */
        static std::vector<std::pair<GroupDim, Int>> coherent_merge(const std::vector<std::pair<GroupDim, Int>> &A,
                                                                    const std::vector<std::pair<GroupDim, Int>> &B) {
            return A.empty() ? B : A;
        }

        /**
         * @brief Projects `AP` onto the geometric entity identified by (`ADim`, `ATag`), mirroring
         * `Blocking::classify()`'s own (private) projection logic — reimplemented here since this
         * class only has access to `Blocking`'s public API.
         * @param AModel The geometric model to look the entity up in.
         * @param ADim Dimension of the entity to project onto.
         * @param ATag `entity_tag()` of the entity to project onto.
         * @param AP The point to project.
         * @return The projected point, or `AP` unchanged if no such entity exists.
         */
        static Point3d project_onto_target(const TGeomModel &AModel, GroupDim ADim, Int ATag, const Point3d &AP) {
            Point3d p = AP;
            if (ADim == GroupDim::Dim0) {
                if (const auto *v = AModel.vertex_by_tag(ATag)) v->project(p);
            } else if (ADim == GroupDim::Dim1) {
                if (const auto *c = AModel.curve_by_tag(ATag)) c->project(p);
            } else if (ADim == GroupDim::Dim2) {
                if (const auto *s = AModel.surface_by_tag(ATag)) s->project(p);
            } else if (ADim == GroupDim::Dim3) {
                if (const auto *vol = AModel.volume_by_tag(ATag)) vol->project(p);
            }
            return p;
        }

        /**
         * @brief Collapses one degenerate 2-edge (bigon) face: merges its 2 duplicate edges into
         * one, removes the bigon, and (if both edges had an "outside" neighbor face) sews those 2
         * neighbors together across the merged edge.
         *
         * A bigon's 2 edges may each have 0 or 1 "outside" neighbor (a face other than the bigon
         * itself). If neither does, both edges — and, transitively, the bigon's whole side of the
         * structure — are garbage-collected by CGAL's own attribute reference-counting once the
         * bigon is removed (same mechanism `Blocking::delete_face()` already relies on); nothing
         * survives, matching `delete_face()`'s own "standalone" precedent. If exactly one does, it
         * survives untouched (its only face is now that one outside neighbor). If both do, `sew<2>`
         * connects the 2 outside neighbors directly — which, per this module's existing sewing
         * behavior (see `Blocking::sew_free_quad_edges()`), also merges the 2 edge attributes into
         * one as a side effect.
         *
         * @param ABigonEntryDart A dart of the bigon face (one of its 2 edges' bigon-side dart).
         * @param AFacesToRebuild Appended with every "outside" neighbor face whose boundary now
         *        includes the surviving merged edge, for the caller to `rebuild_face_geometry()`.
         * @return The surviving edge, or `cmap().null_descriptor` if nothing survived.
         */
        Edge collapse_bigon(Dart ABigonEntryDart, std::vector<Face> &AFacesToRebuild) {
            auto &cmap = m_blocking.cmap();
            const Dart da = ABigonEntryDart;
            const Dart db = cmap.template beta<1>(da);

            const bool a_has_outside = !cmap.template is_free<2>(da);
            const bool b_has_outside = !cmap.template is_free<2>(db);
            const Dart a_outside = a_has_outside ? cmap.template beta<2>(da) : cmap.null_dart_descriptor;
            const Dart b_outside = b_has_outside ? cmap.template beta<2>(db) : cmap.null_dart_descriptor;

            cmap.template remove_cell<2>(da);

            if (a_has_outside && b_has_outside) {
                cmap.template sew<2>(a_outside, b_outside);
                AFacesToRebuild.push_back(cmap.template attribute<2>(a_outside));
                AFacesToRebuild.push_back(cmap.template attribute<2>(b_outside));
                return cmap.template attribute<1>(a_outside);
            }
            if (a_has_outside) {
                AFacesToRebuild.push_back(cmap.template attribute<2>(a_outside));
                return cmap.template attribute<1>(a_outside);
            }
            if (b_has_outside) {
                AFacesToRebuild.push_back(cmap.template attribute<2>(b_outside));
                return cmap.template attribute<1>(b_outside);
            }
            return cmap.null_descriptor;
        }

        /**
         * @brief Rebuilds a surviving merged edge's curve as a straight line between its current
         * 2 endpoints (see the class doc comment for why the 2 pre-merge curves aren't blended).
         * @param AEdge The edge to rebuild.
         * @param ACmap The map `AEdge` belongs to.
         * @return The rebuilt straight curve.
         */
        static TEdgeCurve rebuilt_curve(Edge AEdge, Map &ACmap) {
            const Dart d = AEdge->dart();
            const Point3d &p0 = ACmap.template attribute<0>(d)->info().point;
            const Point3d &p1 = ACmap.template attribute<0>(ACmap.template beta<1>(d))->info().point;
            return BlockingT::straight_curve(p0, p1);
        }

        /** @brief The blocking this operation acts on. */
        BlockingT &m_blocking;
    };

} // namespace gecko
