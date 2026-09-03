#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <vector>

#include <gecko/block/Blocking.h>

namespace gecko {

    /**
     * @class Smoother
     * @brief Smart Laplacian smoothing of a `Blocking`: every corner is pulled towards the average
     * of its neighbours, kept on whatever it is classified on, and actually moved only where that
     * makes the shapes around it better.
     *
     * An operation class rather than a `Blocking` member, for the reason the other multi-step
     * operations are: it carries state across its steps — which corners are locked, who neighbours
     * whom, which cells each corner is measured by — and a `Blocking` is a shape, not a record of
     * what a caller means to do to it. Built against one blocking, kept for as long as the caller
     * wants to go on smoothing it, and holding nothing that outlives an edit to it (see
     * `smooth()`'s note on when to rebuild one).
     *
     * What it does *not* touch is as much the point as what it does:
     *
     * - **Only corner positions.** Control points are never smoothed. They are re-derived once at
     *   the end of a `smooth()` call, from the corners' final positions and each cell's existing
     *   classification — a curved edge that had been fitted onto a model curve is fitted onto it
     *   again, not straightened.
     * - **Never classification.** A corner on a model vertex does not move at all; one on a curve
     *   stays on that curve, one on a surface on that surface. Re-deciding what a cell lies on is
     *   `Blocking::classify()`'s job, and smoothing a classified blocking must not quietly redo it.
     * - **Never the topology.** No cell is created, split or removed, so ids and handles taken
     *   before a `smooth()` call are all still good after it.
     *
     * It works in 2 passes, and the second exists because the first has a ceiling. A smart
     * Laplacian only ever offers a corner *one* position — the average of its neighbours — and takes
     * it or leaves it. Where that one position happens to trade one corner's quality against
     * another's, it is refused, and the blocking settles better than it started and short of
     * regular with every corner refusing to be the one to go first. Nothing about that is a bad
     * average: the right position was simply somewhere between the corner's own and the average,
     * and no rule for accepting or rejecting a single candidate can find it. So the second pass
     * stops proposing and starts *searching* — see `optimize_node()`. `Strategy` picks which passes
     * run; both do by default, the Laplacian first because it is far cheaper and does most of the
     * work in one sweep.
     *
     * @tparam TGeomModel The geometric model type the blocking is built against.
     * @tparam TEdgeCurve The blocking's edge curve representation.
     */
    template<GeomModelConcept TGeomModel, EdgeCurveConcept TEdgeCurve = BezierCurve<Point3d>>
    class Smoother {
    public:
        /** @brief The blocking type this smoother works on. */
        using BlockingT = Blocking<TGeomModel, TEdgeCurve>;
        /** @brief Handle to a node (0-cell) attribute of that blocking. */
        using Node = typename BlockingT::Node;
        /** @brief Handle to a face (2-cell) attribute of that blocking. */
        using Face = typename BlockingT::Face;
        /** @brief Handle to a block (3-cell) attribute of that blocking. */
        using Block = typename BlockingT::Block;

        /** @brief Which of the 2 passes a `smooth()` call runs. */
        enum class Strategy {
            Laplacian,    ///< The smart Laplacian alone — one candidate per corner, taken or left.
            Optimization, ///< The max-min direct search alone, on the blocking as it stands.
            Both          ///< Laplacian first, then the search on what it left behind. The default.
        };

        /** @brief What one `smooth()` call did. */
        struct Report {
            /** @brief Laplacian passes run — fewer than asked for when the blocking settled. */
            std::size_t laplacian_passes = 0;
            /** @brief Optimization passes run, likewise. */
            std::size_t optimization_passes = 0;
            /** @brief Corners moved, summed over every pass (a corner moved twice counts twice). */
            std::size_t moves = 0;
            /** @brief The worst cell quality in the blocking afterwards, as `Blocking::block_quality()`
             * and `Blocking::face_quality()` measure it. Never lower than it was before the call. */
            double worst_quality = 0.0;
        };

        /**
         * @brief Builds a smoother for one blocking.
         * @param ABlocking The blocking to smooth; must outlive this smoother (only a non-owning
         *        pointer is stored).
         */
        explicit Smoother(BlockingT &ABlocking) : m_blocking(&ABlocking) {}

        /**
         * @brief Locks one corner, so no `smooth()` call moves it.
         *
         * By id and not by handle, because a lock is meant to outlive the edits made between 2
         * smoothing passes and a handle is not: CGAL discards and rebuilds an attribute whenever the
         * orbit behind it is disturbed, which the ordinary editing operations do all the time (see
         * `CellInfo::id`). An id that no longer names a node is simply never matched.
         *
         * @param ANodeId The node's `CellInfo::id`.
         */
        void lock(Int ANodeId) { m_locked.insert(ANodeId); }

        /** @brief Unlocks one corner. @param ANodeId The node's `CellInfo::id`. */
        void unlock(Int ANodeId) { m_locked.erase(ANodeId); }

        /** @brief Replaces the whole set of locked corners. @param ANodeIds Their `CellInfo::id`s. */
        void set_locked(const std::vector<Int> &ANodeIds) {
            m_locked = std::set<Int>(ANodeIds.begin(), ANodeIds.end());
        }

        /** @brief Whether a corner is locked. @param ANodeId The node's `CellInfo::id`. @return true if it is. */
        [[nodiscard]] bool is_locked(Int ANodeId) const { return m_locked.contains(ANodeId); }

        /** @brief The locked corners. @return Their `CellInfo::id`s. */
        [[nodiscard]] const std::set<Int> &locked() const { return m_locked; }

        /**
         * @brief Smooths the blocking, at most @p AIterations times over.
         *
         * One iteration walks every corner in id order and, for each, averages the neighbours that
         * are *entitled to have a say* (see `laplacian_target()`), pulls that average back onto
         * whatever the corner is classified on, and keeps it only if the worst cell around the
         * corner is no worse than it was. That last test is what makes this the "smart" Laplacian
         * rather than the plain one, and it is what stops the classic failure: plain Laplacian
         * smoothing happily drags a corner across its own neighbours and turns a block inside out,
         * which no amount of iterating undoes. Since every accepted move leaves the cells around it
         * no worse and touches no others, the blocking's own worst cell can only improve across a
         * whole call.
         *
         * Corner by corner and in id order, not all at once from a snapshot. Taking the whole
         * iteration from the positions it started at would make the result independent of the order
         * corners are visited in, which is tempting — but it also breaks the guarantee above, since
         * a move measured as an improvement against the old positions need not be one once every
         * other corner has moved too. Walking in id order gets determinism the other way: ids are
         * stable names, so the same blocking smooths to the same result, and each move is measured
         * against the positions it is actually made against.
         *
         * Iterating stops early the moment a whole pass moves nothing — from there no later pass
         * would either, every corner having been offered its move and refused it.
         *
         * The optimization passes then run over the same corners in the same order, each one
         * searching for the position that leaves the worst cell around it as good as possible
         * (`optimize_node()`), and stopping early on the same rule. @p AIterations caps each kind of
         * pass separately: a call is allowed that many Laplacian passes and that many optimization
         * passes, not that many between them.
         *
         * Control points are re-derived once, at the end, and only around the corners that actually
         * moved (`Blocking::refit_around()`). Not per iteration, which would rebuild the same cells
         * over and over to feed positions nothing reads in between: the quality test measures
         * corners, deliberately, so that intermediate geometry is never needed.
         *
         * @param AIterations The maximum number of passes of each kind.
         * @param AStrategy Which passes to run; both by default.
         * @return What the call did, and the blocking's worst cell quality afterwards.
         * @pre The blocking's nodes must have been named (`CellInfo::id` >= 0), which every way of
         *      creating them does.
         * @note Holds no state across calls beyond the locked set, so an edit to the blocking
         *       between 2 calls is picked up: adjacency and incidence are rebuilt by each call.
         */
        Report smooth(std::size_t AIterations, Strategy AStrategy = Strategy::Both) {
            Report report;

            const std::vector<Node> order = nodes_in_id_order();
            const std::map<Node, std::vector<Node>> neighbours = build_adjacency();
            const std::map<Node, std::vector<Block>> blocks = build_block_incidence();
            const std::map<Node, std::vector<Face>> faces = build_standalone_face_incidence();

            std::set<Node> moved;
            if (AStrategy != Strategy::Optimization) {
                for (std::size_t pass = 0; pass < AIterations; ++pass) {
                    std::size_t moves = 0;
                    for (const Node node : order) {
                        if (try_move(node, neighbours, blocks, faces)) {
                            moved.insert(node);
                            ++moves;
                        }
                    }
                    report.laplacian_passes = pass + 1;
                    report.moves += moves;
                    if (moves == 0) break;
                }
            }
            if (AStrategy != Strategy::Laplacian) {
                for (std::size_t pass = 0; pass < AIterations; ++pass) {
                    std::size_t moves = 0;
                    for (const Node node : order) {
                        if (optimize_node(node, neighbours, blocks, faces)) {
                            moved.insert(node);
                            ++moves;
                        }
                    }
                    report.optimization_passes = pass + 1;
                    report.moves += moves;
                    if (moves == 0) break;
                }
            }

            if (!moved.empty()) {
                std::vector<Point3d> points;
                points.reserve(moved.size());
                for (const Node node : moved) {
                    points.push_back(node->info().point);
                }
                // Infer-only, so the tolerances are never read: every cell here takes the
                // classification its own boundary agrees on, and a cell whose boundary agrees on
                // nothing is left unclassified rather than searched for something nearby.
                m_blocking->refit_around(points, typename BlockingT::Tolerances{});
            }

            report.worst_quality = worst_quality();
            return report;
        }

        /**
         * @brief The worst cell quality anywhere in the blocking — its hex blocks measured by
         * `Blocking::block_quality()`, its standalone quads by `Blocking::face_quality()`.
         * @return That minimum, or 1 for a blocking with no cell to measure.
         */
        double worst_quality() {
            double worst = std::numeric_limits<double>::max();
            auto &map = m_blocking->cmap();
            for (auto it = map.template attributes<3>().begin(), end = map.template attributes<3>().end(); it != end;
                 ++it) {
                const Block block = it;
                worst = std::min(worst, m_blocking->block_quality(block));
            }
            for (auto it = map.template attributes<2>().begin(), end = map.template attributes<2>().end(); it != end;
                 ++it) {
                const Face face = it;
                if (!m_blocking->belongs_to_block(face)) worst = std::min(worst, m_blocking->face_quality(face));
            }
            return (worst == std::numeric_limits<double>::max()) ? 1.0 : worst;
        }

    private:
        /**
         * @brief Offers one corner its Laplacian move and takes it only if it is an improvement.
         * @param ANode The corner.
         * @param ANeighbours The whole adjacency table.
         * @param ABlocks Which blocks each corner is measured by.
         * @param AFaces Which standalone faces each corner is measured by.
         * @return true if the corner moved.
         */
        bool try_move(Node ANode,
                      const std::map<Node, std::vector<Node>> &ANeighbours,
                      const std::map<Node, std::vector<Block>> &ABlocks,
                      const std::map<Node, std::vector<Face>> &AFaces) {
            if (is_locked(ANode->info().id)) return false;
            if (pinned_to_vertex(ANode)) return false;

            const auto found = ANeighbours.find(ANode);
            if (found == ANeighbours.end()) return false;

            const auto target = laplacian_target(ANode, found->second);
            if (!target.has_value()) return false;

            const Point3d before = ANode->info().point;
            const Point3d after = m_blocking->project_onto_classification(ANode, *target);
            if (!worth_moving(before, after, found->second)) return false;

            const double was = local_quality(ANode, ABlocks, AFaces);
            ANode->info().point = after;
            const double now = local_quality(ANode, ABlocks, AFaces);
            // Kept unless it makes something worse — not "kept only if it makes something better".
            // The 2 are not the same rule and the strict one does not work: the corner fixing the
            // local minimum is usually not the one being offered the move, so a move that tidies up
            // everything else while leaving that corner exactly as it was reads as no improvement
            // and is refused. A perturbed grid then stalls a long way from regular, each corner
            // waiting for another to go first. Refusing only what actually degrades keeps the whole
            // point of the test — nothing is ever made worse, and in particular nothing is ever
            // turned inside out — while letting the pass make progress.
            if (now >= was) return true;

            ANode->info().point = before;
            return false;
        }

        /**
         * @brief Searches for the position that leaves the worst cell around one corner as good as
         * it can be, and moves the corner there.
         *
         * Where `try_move()` proposes a single position and takes it or leaves it, this looks. It is
         * a compass search: from the corner's current position, try a step along each of the 6 axis
         * directions, keep any that leaves the worst cell around the corner better than the best
         * found so far, and when a whole round of 6 finds nothing, halve the step and go round
         * again — down to a step small against the corner's own neighbourhood. Every candidate is
         * pulled onto the corner's classification before being measured, so the search runs *on*
         * the curve or surface the corner belongs to rather than in space with a correction
         * afterwards, and the axis directions need no relation to that entity's own parameters.
         *
         * Direct search and not a gradient method, deliberately: the thing being maximized is the
         * *minimum* over the incident cells, which has a corner exactly where the cell attaining
         * that minimum changes hands — precisely where the interesting positions are. A gradient
         * there is either undefined or, worse, defined and misleading. Trying directions asks no
         * question the objective cannot answer.
         *
         * The corner's current position is in the set being searched, so the answer is never worse
         * than doing nothing: this pass, like the Laplacian one, can only raise the worst cell in
         * the blocking.
         *
         * @param ANode The corner.
         * @param ANeighbours The whole adjacency table, which sets the search's length scale.
         * @param ABlocks Which blocks each corner is measured by.
         * @param AFaces Which standalone faces each corner is measured by.
         * @return true if the corner moved.
         */
        bool optimize_node(Node ANode,
                           const std::map<Node, std::vector<Node>> &ANeighbours,
                           const std::map<Node, std::vector<Block>> &ABlocks,
                           const std::map<Node, std::vector<Face>> &AFaces) {
            if (is_locked(ANode->info().id)) return false;
            if (pinned_to_vertex(ANode)) return false;

            const auto found = ANeighbours.find(ANode);
            if (found == ANeighbours.end()) return false;

            const Point3d origin = ANode->info().point;
            double scale = 0.0;
            for (const Node neighbour : found->second) {
                scale = std::max(scale, Vector3d(origin, neighbour->info().point).norm());
            }
            if (scale == 0.0) return false;

            double best_quality = local_quality(ANode, ABlocks, AFaces);
            if (best_quality == std::numeric_limits<double>::max()) return false;

            Point3d best = origin;
            double step = 0.25 * scale;
            const double smallest = 1e-4 * scale;
            // Improvements have to clear a margin, and the rounds are capped. Both guard the same
            // thing: a search that goes on trading gains of 1e-17 for evaluations forever, which is
            // what a strict `>` on a quantity computed from square roots eventually finds.
            const double margin = 1e-12;
            for (std::size_t round = 0; round < MAX_SEARCH_ROUNDS && step > smallest; ++round) {
                bool improved = false;
                for (const Vector3d &direction : SEARCH_DIRECTIONS) {
                    const Point3d candidate = m_blocking->project_onto_classification(ANode, best + direction * step);
                    ANode->info().point = candidate;
                    const double quality = local_quality(ANode, ABlocks, AFaces);
                    if (quality > best_quality + margin) {
                        best_quality = quality;
                        best = candidate;
                        improved = true;
                    }
                }
                if (!improved) step *= 0.5;
            }

            ANode->info().point = best;
            return best != origin;
        }

        /**
         * @brief Where the Laplacian would put one corner: the average of the neighbours entitled to
         * have a say in it.
         *
         * Which neighbours those are is what makes this classification-aware. A corner lying on a
         * model curve is averaged over the neighbours that lie on that same curve — or on one of its
         * end vertices — and not over the ones heading off into the interior, which would drag it
         * along the curve towards the middle of the blocking before the projection pulled it back
         * on. Same one dimension up: a corner on a surface hears from the neighbours on that surface
         * and on the curves bounding it, and not from the interior.
         *
         * Stated once and for all dimensions at once through the model's own containment: a
         * neighbour has a say when the corner's own classification target is among the entities
         * *containing* the neighbour's (`Blocking::containing_set()`, which counts an entity as
         * containing itself). A corner in the interior — classified on a volume, or not classified
         * at all — hears from all of its neighbours, there being nothing to stay on.
         *
         * @param ANode The corner.
         * @param ANeighbours Its neighbours.
         * @return Their average, or nothing if none of them has a say.
         */
        std::optional<Point3d> laplacian_target(Node ANode, const std::vector<Node> &ANeighbours) const {
            const auto &targets = ANode->info().geom_targets;
            const bool interior = targets.empty() || lowest_dim(targets) == GroupDim::Dim3;

            Vector3d sum;
            std::size_t count = 0;
            for (const Node neighbour : ANeighbours) {
                if (!interior) {
                    const auto containing = m_blocking->containing_set(neighbour->info().geom_targets);
                    if (std::ranges::none_of(targets,
                                             [&containing](const auto &t) { return containing.contains(t); })) {
                        continue;
                    }
                }
                sum += Vector3d(ANode->info().point, neighbour->info().point);
                ++count;
            }
            if (count == 0) return std::nullopt;
            return ANode->info().point + sum / static_cast<double>(count);
        }

        /**
         * @brief Whether a proposed move is big enough to be worth making.
         *
         * Measured against how far the corner's own neighbours are, not against a fixed distance: a
         * blocking is whatever size its model is, and an absolute threshold would be everything on
         * one and nothing on another. Below it the move is refused outright, which keeps a corner
         * already sitting where the Laplacian wants it from being offered — and charged for — a
         * rounding-error move on every pass forever.
         *
         * @param ABefore Where the corner is.
         * @param AAfter Where it would go.
         * @param ANeighbours Its neighbours, which set the scale.
         * @return true if the move is worth making.
         */
        static bool worth_moving(const Point3d &ABefore, const Point3d &AAfter, const std::vector<Node> &ANeighbours) {
            double scale = 0.0;
            for (const Node neighbour : ANeighbours) {
                scale = std::max(scale, Vector3d(ABefore, neighbour->info().point).norm());
            }
            return Vector3d(ABefore, AAfter).norm() > 1e-12 * std::max(scale, 1.0);
        }

        /**
         * @brief The worst quality among the cells one corner belongs to.
         * @param ANode The corner.
         * @param ABlocks Which blocks each corner is measured by.
         * @param AFaces Which standalone faces each corner is measured by.
         * @return That minimum, or +infinity for a corner belonging to no cell at all — which no
         *         move can improve on, so such a corner stays put.
         */
        double local_quality(Node ANode,
                             const std::map<Node, std::vector<Block>> &ABlocks,
                             const std::map<Node, std::vector<Face>> &AFaces) {
            double worst = std::numeric_limits<double>::max();
            if (const auto found = ABlocks.find(ANode); found != ABlocks.end()) {
                for (const Block block : found->second) {
                    worst = std::min(worst, m_blocking->block_quality(block));
                }
            }
            if (const auto found = AFaces.find(ANode); found != AFaces.end()) {
                for (const Face face : found->second) {
                    worst = std::min(worst, m_blocking->face_quality(face));
                }
            }
            return worst;
        }

        /** @brief Whether a corner sits on a model vertex, and so cannot move at all.
         * @param ANode The corner. @return true if it is classified on a vertex. */
        static bool pinned_to_vertex(Node ANode) {
            const auto &targets = ANode->info().geom_targets;
            return !targets.empty() && lowest_dim(targets) == GroupDim::Dim0;
        }

        /** @brief The lowest dimension a cell is classified on. @param ATargets Its `geom_targets`,
         * which must not be empty. @return That dimension. */
        static GroupDim lowest_dim(const std::vector<std::pair<GroupDim, Int>> &ATargets) {
            GroupDim lowest = GroupDim::Undefined;
            for (const auto &[dim, tag] : ATargets) {
                if (dim < lowest) lowest = dim;
            }
            return lowest;
        }

        /**
         * @brief Every corner of the blocking, in `CellInfo::id` order.
         * @return The corners, ordered — which is what makes a smoothing pass reproducible, the
         *         order `attributes<0>()` traverses in being CGAL's to change.
         */
        std::vector<Node> nodes_in_id_order() const {
            std::vector<Node> nodes;
            auto &map = m_blocking->cmap();
            for (auto it = map.template attributes<0>().begin(), end = map.template attributes<0>().end(); it != end;
                 ++it) {
                nodes.push_back(it);
            }
            std::ranges::sort(nodes, [](Node AA, Node AB) { return AA->info().id < AB->info().id; });
            return nodes;
        }

        /**
         * @brief Who neighbours whom, along the blocking's edges.
         *
         * Read off the edge attributes rather than walked from each corner's own dart orbit, for the
         * reason `Blocking::move_node()` documents: on a blocking whose blocks are not all sewn
         * together, a corner's orbit is a single dart and the incidence iterators miss most of what
         * actually touches it.
         *
         * @return One entry per corner that has at least one neighbour.
         */
        std::map<Node, std::vector<Node>> build_adjacency() const {
            std::map<Node, std::vector<Node>> adjacency;
            auto &map = m_blocking->cmap();
            for (auto it = map.template attributes<1>().begin(), end = map.template attributes<1>().end(); it != end;
                 ++it) {
                const auto dart = it->dart();
                const Node n0 = map.template attribute<0>(dart);
                const Node n1 = map.template attribute<0>(map.template beta<1>(dart));
                if (n0 == n1) continue;
                adjacency[n0].push_back(n1);
                adjacency[n1].push_back(n0);
            }
            // Deduplicated: 2 corners can be joined by more than one edge, and a neighbour counted
            // twice would weigh twice as much in an average that is meant to be over neighbours.
            for (auto &[node, list] : adjacency) {
                std::ranges::sort(list);
                list.erase(std::ranges::unique(list).begin(), list.end());
            }
            return adjacency;
        }

        /** @brief Which hex blocks each corner belongs to. @return One entry per corner of a block. */
        std::map<Node, std::vector<Block>> build_block_incidence() const {
            std::map<Node, std::vector<Block>> incidence;
            auto &map = m_blocking->cmap();
            for (auto it = map.template attributes<3>().begin(), end = map.template attributes<3>().end(); it != end;
                 ++it) {
                const Block block = it;
                for (const Node corner : m_blocking->frame_of(block)) {
                    incidence[corner].push_back(block);
                }
            }
            return incidence;
        }

        /** @brief Which standalone (2D block) faces each corner belongs to — a face bounding a hex
         * is left out, its block measuring it better. @return One entry per corner of such a face. */
        std::map<Node, std::vector<Face>> build_standalone_face_incidence() const {
            std::map<Node, std::vector<Face>> incidence;
            auto &map = m_blocking->cmap();
            for (auto it = map.template attributes<2>().begin(), end = map.template attributes<2>().end(); it != end;
                 ++it) {
                const Face face = it;
                if (m_blocking->belongs_to_block(face)) continue;
                for (const Node corner : m_blocking->frame_of(face)) {
                    incidence[corner].push_back(face);
                }
            }
            return incidence;
        }

        /** @brief The directions `optimize_node()` tries, one step along each per round: the 6 axis
         * directions. Axis-aligned and not aligned to anything the corner lies on, because every
         * candidate is projected onto its classification anyway — a step normal to a curve simply
         * comes back where it started and costs one evaluation. */
        static constexpr std::array<Vector3d, 6> SEARCH_DIRECTIONS = {Vector3d(1.0, 0.0, 0.0),
                                                                      Vector3d(-1.0, 0.0, 0.0),
                                                                      Vector3d(0.0, 1.0, 0.0),
                                                                      Vector3d(0.0, -1.0, 0.0),
                                                                      Vector3d(0.0, 0.0, 1.0),
                                                                      Vector3d(0.0, 0.0, -1.0)};

        /** @brief How many rounds of `optimize_node()`'s compass search one corner gets, however the
         * step size is going. A backstop, not a tuning knob: the step-size floor ends the search
         * long before this on anything well behaved. */
        static constexpr std::size_t MAX_SEARCH_ROUNDS = 200;

        /** @brief The blocking being smoothed, which this smoother does not own. */
        BlockingT *m_blocking;
        /** @brief The corners no pass may move, by `CellInfo::id`. */
        std::set<Int> m_locked;
    };

} // namespace gecko
