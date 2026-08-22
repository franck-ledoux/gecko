#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <vector>

#include <gecko/block/CellData.h>
#include <gecko/geom_itf/Concepts.h>
#include <gecko/math/CoonsPatch.h>
#include <gecko/math/CurveSurfaceTraits.h>
#include <gecko/math/Point3d.h>
#include <gecko/math/Vector3d.h>
#include <gecko/mesh/UnstructuredMesh.h>

namespace gecko {

    /**
     * @class Blocking
     * @brief A structured (quad/hex) blocking of a geometric model, backed by a single,
     * always-dimension-3 `CGAL::Combinatorial_map`.
     *
     * A "2D" (quad) block is a standalone 2-cell with no incident 3-cell; a "3D" (hex) block is a
     * genuine 3-cell — both live in, and can be freely mixed within, the same map (see CellData.h).
     * `TEdgeCurve` (default `BezierCurve<Point3d>`) names the *representation family* only — the
     * degree is a property of the blocking, set at construction and changeable afterwards through
     * `set_degree()`. Any degree gives a curved blocking and degree 1 collapses exactly to linear
     * blocking, through the same Coons/TFI construction (see math/CoonsPatch.h). Keeping the family
     * a template parameter is what leaves room for a future `NurbsCurve` to take its place without
     * touching this class; keeping the degree out of it is what lets a structure being edited change
     * order in place, which a degree baked into the type could not do without rebuilding every cell
     * into a different C++ type. `TGeomModel` (`GeomModelConcept`-conforming, e.g. `FacetedGeometry`)
     * is the geometric model this blocking is built against, held for the whole lifetime of the
     * object — "a block structure always exists along a geometric model".
     *
     * @tparam TGeomModel Geometric model type, must satisfy `GeomModelConcept`.
     * @tparam TEdgeCurve Edge curve representation, must satisfy `EdgeCurveConcept`.
     */
    template<GeomModelConcept TGeomModel, EdgeCurveConcept TEdgeCurve = BezierCurve<Point3d>>
    class Blocking {
        /** @brief The face representation paired with `TEdgeCurve` via `CurveSurfaceTraits` (same
         * order per issue #22) — `BezierSurface<N,Point3d>` for a `BezierCurve<N,Point3d>` edge
         * today; a future `NurbsCurve`/`BSplineCurve` brings its own paired surface with it, with no
         * change needed here. */
        using FaceSurfaceT = typename CurveSurfaceTraits<TEdgeCurve>::Surface;
        /** @brief The block representation paired with `FaceSurfaceT` via `SurfaceVolumeTraits`
         * (same order per issue #22) — `BezierVolume<N,Point3d>` for a `BezierSurface<N,Point3d>`
         * face today. */
        using BlockVolumeT = typename SurfaceVolumeTraits<FaceSurfaceT>::Volume;

    public:
        /** @brief The underlying combinatorial map type (always dimension 3). */
        using Map = CMap<TEdgeCurve, FaceSurfaceT, BlockVolumeT>;
        /** @brief A dart descriptor into the underlying map. */
        using Dart = typename Map::Dart_descriptor;
        /** @brief Handle to a block (3-cell) attribute. */
        using Block = typename Map::template Attribute_handle<3>::type;
        /** @brief Handle to a face (2-cell) attribute. */
        using Face = typename Map::template Attribute_handle<2>::type;
        /** @brief Handle to an edge (1-cell) attribute. */
        using Edge = typename Map::template Attribute_handle<1>::type;
        /** @brief Handle to a node (0-cell) attribute. */
        using Node = typename Map::template Attribute_handle<0>::type;

        /**
         * @brief Name of the per-node `Variable<Int>` `to_mesh()` writes giving the dimension of the
         * geometric entity each mesh node is classified on: 0 vertex, 1 curve, 2 surface, 3 volume,
         * or -1 if unclassified. Every mesh node inherits the classification of the block structure
         * cell it was sampled from — a corner node's own, an edge/face/block-interior point's
         * generating edge/face/block (see `CellData.h`'s `geom_targets`, and `node_classification_dims()`'s
         * matching "first target speaks for all" convention for a node classified on several entities
         * at once).
         */
        static constexpr std::string_view NODE_CLASSIFICATION_DIM_VARIABLE = "classification_dim";
        /**
         * @brief Name of the per-node `Variable<Int>` `to_mesh()` writes giving the `entity_tag()` of
         * the geometric entity named by `NODE_CLASSIFICATION_DIM_VARIABLE`, or -1 if unclassified.
         */
        static constexpr std::string_view NODE_CLASSIFICATION_TAG_VARIABLE = "classification_tag";

        /**
         * @brief Per-dimension distance thresholds for snapping a blocking onto its geometric model.
         *
         * Three separate values rather than one because the scales genuinely differ: 2 distinct
         * vertices are typically far closer to each other than to any curve, so a tolerance loose
         * enough to catch a surface would snap a corner to the wrong vertex.
         */
        struct Tolerances {
            /** @brief Threshold for snapping onto a vertex (dimension 0). */
            double vertex = 0.0;
            /** @brief Threshold for snapping onto a curve (dimension 1). */
            double curve = 0.0;
            /** @brief Threshold for snapping onto a surface (dimension 2), and onto a volume. */
            double surface = 0.0;
        };

        /**
         * @brief Constructor.
         * @param AGeomModel The geometric model this blocking is built against; must outlive this
         *        Blocking (only a non-owning pointer is stored).
         * @param ADegree The degree every edge curve, face surface and block volume is built at.
         *        1 gives straight edges. Changeable later through `set_degree()`.
         */
        explicit Blocking(const TGeomModel &AGeomModel, std::size_t ADegree = 1)
            : m_degree(ADegree), m_geom_model(&AGeomModel) {
            assert(ADegree >= 1 && "Blocking: degree must be at least 1");
        }

        /** @brief The degree every cell's geometry is currently built at. @return The degree. */
        [[nodiscard]] std::size_t degree() const { return m_degree; }

        /**
         * @brief Rebuilds the whole blocking's geometry at a new degree, refitting it onto the
         * geometric model as it goes.
         *
         * Topology is untouched: the same nodes, edges, faces and blocks come out, carrying the same
         * classification. Only the representation changes — every edge is rebuilt at @p ADegree and
         * refitted onto whatever it is classified on, and every face and block is rebuilt from those
         * edges. Raising the order therefore does not merely add control points, it *uses* them: an
         * edge lying on a curved model curve, which at degree 1 could only be its chord, comes back
         * following it.
         *
         * Implemented as a reclassification rather than as degree elevation, deliberately. Elevation
         * is exact and would preserve each edge's current shape, but that shape is the best fit at
         * the *old* degree — which is precisely what the caller is asking to improve. An edge with
         * nothing to fit onto still ends up where elevation would have put it, since a straight edge
         * rebuilt straight at any degree is the same curve.
         *
         * @param ADegree The new degree, at least 1.
         * @param ATolVertex Tolerance for snapping onto a vertex, as `classify()` takes it.
         * @param ATolCurve Tolerance for snapping onto a curve. Defaults to @p ATolVertex.
         * @param ATolSurface Tolerance for snapping onto a surface. Defaults to the resolved curve
         *        tolerance.
         * @return false, changing nothing, if @p ADegree is less than 1.
         */
        bool set_degree(std::size_t ADegree, double ATolVertex, double ATolCurve = -1.0, double ATolSurface = -1.0) {
            if (ADegree < 1) return false;
            m_degree = ADegree;
            // classify() rebuilds every cell's geometry from scratch — edges refitted onto their own
            // classification, faces and blocks re-derived from those edges — and every one of those
            // rebuilds now happens at the degree just set.
            classify(ATolVertex, ATolCurve, ATolSurface);
            return true;
        }

        /**
         * @brief Creates a new, unsewn, straight-edged quad block (a standalone face, no 3-cell),
         * with its Coons surface built from the 4 straight boundary edges.
         * @param ACorners The 4 corner positions, in perimeter (CCW as seen from outside) order —
         *        matching `CubicTraits::Face::EdgeNodes`' `{0,1},{1,2},{2,3},{3,0}` perimeter.
         * @return The newly created face.
         */
        Face create_quad_block(const std::array<Point3d, 4> &ACorners) {
            std::array<Node, 4> nodes{};
            for (std::size_t i = 0; i < 4; ++i) {
                nodes[i] = create_node(ACorners[i]);
            }
            Dart d1 = m_cmap.make_quadrangle(nodes[0], nodes[1], nodes[2], nodes[3]);

            // The 4 boundary edges, built once, directly in the orientation coons_surface_from_edges
            // needs: EdgeU0 (corner0->corner1), EdgeU1 (corner3->corner2), EdgeV0 (corner0->corner3),
            // EdgeV1 (corner1->corner2) — u along corner0->corner1, v along corner0->corner3.
            std::array<TEdgeCurve, 4> curves{};
            for (std::size_t k = 0; k < 4; ++k) {
                curves[k] = straight_curve(ACorners[static_cast<std::size_t>(QUAD_EDGES[k].first)],
                                           ACorners[static_cast<std::size_t>(QUAD_EDGES[k].second)]);
            }

            for (auto it = m_cmap.template one_dart_per_incident_cell<1, 2>(d1).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<1, 2>(d1).end();
                 it != itend;
                 ++it) {
                const int a = node_index(nodes, m_cmap.template attribute<0>(it));
                const int b = node_index(nodes, m_cmap.template attribute<0>(m_cmap.template beta<1>(it)));
                const std::size_t k = find_edge(QUAD_EDGES, a, b);
                Edge e = create_edge();
                m_cmap.template set_attribute<1>(it, e);
                // `curves[k]` runs from QUAD_EDGES[k].first to .second; store_curve() flips it if
                // this edge's own dart runs the other way, so the class' orientation invariant holds
                // from the moment the edge exists.
                store_curve(e, curves[k], nodes[static_cast<std::size_t>(QUAD_EDGES[k].first)]);
            }

            Face f = m_cmap.template create_attribute<2>();
            m_cmap.template set_attribute<2>(d1, f);
            // Pinned, because `set_attribute` is free to leave the attribute pointing at any dart of
            // the cell, and everything that later reads a face's surface — `build_face_grid()`,
            // `rebuild_face_surface()` — anchors `(u,v)` to the cycle starting at *that* dart. Left
            // to chance, the stored surface and its readers would disagree by a square symmetry.
            f->set_dart(d1);
            // Built from that cycle rather than from `curves` directly, so the frame holds whatever
            // dart the cell ends up on later.
            rebuild_face_surface(f);
            return f;
        }

        /**
         * @brief Creates a new, unsewn, straight-edged hex block (a genuine 3-cell), with its 6
         * faces' Coons surfaces and its own TFI volume built from the 12 straight boundary edges.
         * @param ACorners The 8 corner positions, matching `CubicTraits::Cell`'s own documented
         *        HEX8 node ordering (bottom perimeter 0-1-2-3, top perimeter 4-5-6-7 directly
         *        above) — verified (not assumed) to coincide with this class' own internal
         *        CGAL-hexahedron corner correspondence, see `blocking_creation_tests.cpp`.
         * @return The newly created block.
         */
        Block create_hex_block(const std::array<Point3d, 8> &ACorners) {
            Dart d1 = m_cmap.make_combinatorial_hexahedron();

            std::array<Node, 8> nodes{};
            nodes[0] = create_node(ACorners[0]);
            m_cmap.template set_attribute<0>(d1, nodes[0]);
            nodes[1] = create_node(ACorners[1]);
            m_cmap.template set_attribute<0>(m_cmap.template beta<1>(d1), nodes[1]);
            nodes[2] = create_node(ACorners[2]);
            m_cmap.template set_attribute<0>(m_cmap.template beta<1, 1>(d1), nodes[2]);
            nodes[3] = create_node(ACorners[3]);
            m_cmap.template set_attribute<0>(m_cmap.template beta<1, 1, 1>(d1), nodes[3]);
            nodes[4] = create_node(ACorners[4]);
            m_cmap.template set_attribute<0>(m_cmap.template beta<2, 1, 1>(d1), nodes[4]);
            nodes[5] = create_node(ACorners[5]);
            m_cmap.template set_attribute<0>(m_cmap.template beta<1, 2, 1, 1>(d1), nodes[5]);
            nodes[6] = create_node(ACorners[6]);
            m_cmap.template set_attribute<0>(m_cmap.template beta<1, 1, 2, 1, 1>(d1), nodes[6]);
            nodes[7] = create_node(ACorners[7]);
            m_cmap.template set_attribute<0>(m_cmap.template beta<1, 1, 1, 2, 1, 1>(d1), nodes[7]);

            // The 12 boundary edges, built once from ACorners, in the fixed direction HEX_FACES
            // below always references them in (verified by hand — see blocking_creation_tests.cpp
            // — to be used identically by both of each edge's 2 incident faces, so no per-face
            // reversal is ever needed).
            std::array<TEdgeCurve, 12> curves{};
            for (std::size_t k = 0; k < 12; ++k) {
                curves[k] = straight_curve(ACorners[static_cast<std::size_t>(HEX_EDGES[k].first)],
                                           ACorners[static_cast<std::size_t>(HEX_EDGES[k].second)]);
            }

            for (auto it = m_cmap.template one_dart_per_incident_cell<1, 3>(d1).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<1, 3>(d1).end();
                 it != itend;
                 ++it) {
                const int a = node_index(nodes, m_cmap.template attribute<0>(it));
                const int b = node_index(nodes, m_cmap.template attribute<0>(m_cmap.template beta<1>(it)));
                const std::size_t k = find_edge(HEX_EDGES, a, b);
                Edge e = create_edge();
                m_cmap.template set_attribute<1>(it, e);
                store_curve(e, curves[k], nodes[static_cast<std::size_t>(HEX_EDGES[k].first)]);
            }

            // The 6 bounding faces' Coons surfaces, built once from the 12 edges above, in the
            // Fu0/Fu1/Fv0/Fv1/Fw0/Fw1 order tfi_volume_from_faces() expects.
            std::array<FaceSurfaceT, 6> face_surfaces{};
            for (std::size_t f = 0; f < 6; ++f) {
                const auto &spec = HEX_FACES[f];
                face_surfaces[f] = coons_surface_from_edges(
                    curves[spec.edge_u0], curves[spec.edge_u1], curves[spec.edge_v0], curves[spec.edge_v1]);
            }

            for (auto it = m_cmap.template one_dart_per_incident_cell<2, 3>(d1).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<2, 3>(d1).end();
                 it != itend;
                 ++it) {
                std::array<int, 4> corner_idx{};
                Dart fd = it;
                for (std::size_t c = 0; c < 4; ++c) {
                    corner_idx[c] = node_index(nodes, m_cmap.template attribute<0>(fd));
                    fd = m_cmap.template beta<1>(fd);
                }
                const std::size_t f = find_face(corner_idx);
                (void)f;
                Face face = create_face();
                m_cmap.template set_attribute<2>(it, face);
                face->set_dart(it);
                // Deliberately *not* `face_surfaces[f]`: that one is parameterized in HEX_FACES'
                // frame, while everything that later reads a face's surface — `build_face_grid()`,
                // `rebuild_face_surface()` — anchors `(u,v)` to the face's own dart cycle, which is
                // some rotation or reflection of it. Storing the HEX_FACES form would leave the 2
                // disagreeing by a square symmetry. That went unnoticed for as long as faces were
                // only ever sampled on a uniform grid, which a transpose maps onto itself; a sheet
                // cut is the first thing to sample one asymmetrically, and it put the edge it
                // inserts across the wrong axis. Rebuilding here yields the very same Coons patch,
                // in the frame the rest of the class agrees on.
                rebuild_face_surface(face);
                m_hex_faces.push_back(face);
            }

            Block b = m_cmap.template create_attribute<3>();
            m_cmap.template set_attribute<3>(d1, b);
            // Pinned for the same reason as the faces above: `emit_hex_cells()` and
            // `rebuild_block_volume()` both walk the block's corners from whatever dart the
            // attribute holds, so leaving it unpinned lets the stored volume and its readers
            // disagree by a cube symmetry.
            b->set_dart(d1);
            // Rebuilt from the block's own dart cycle rather than assembled from `face_surfaces`
            // directly, for the same reason the faces above are: the TFI of those 6 surfaces is
            // parameterized in the order the caller passed its corners in, while `emit_hex_cells()`
            // and `rebuild_block_volume()` both anchor `(u,v,w)` to the cycle CGAL gives the 3-cell.
            // The 2 differ by a cube symmetry whenever those disagree, which a uniform sampling
            // hides and a sheet cut does not — it would split the block across the wrong axis.
            rebuild_block_volume(b);
            return b;
        }

        /**
         * @brief Auto-detects coincident corner positions across every block created so far (any
         * mix of 2D quads and 3D hexes) and sews them, gluing hex blocks at shared faces
         * (dimension-3 sew, only ever considered between 2 genuine hex-boundary faces) and quad
         * blocks at shared edges (dimension-2 sew, only ever considered between 2 standalone quad
         * edges — a hex's own 12 edges are already sewn to their 2 sibling faces *within* that same
         * hex at creation time, so they're never free at dimension 2 and never candidates here).
         * Not incremental: call again after adding more blocks to pick up new adjacencies. Mirrors
         * `UnstructuredMesh::build_connectivity()`'s "create raw entities, then link them" pattern.
         */
        void build_connectivity() {
            sew_matching<3>(m_hex_faces, 4);
            sew_free_quad_edges();
        }

        /**
         * @brief Classifies every node/edge/face/block of the blocking onto `geom_model()`'s
         * vertices/curves/surfaces/volumes, then refits/rebuilds their geometry to conform.
         *
         * **Nodes** are classified by proximity: the search runs dimension by dimension and stops at
         * the lowest one with any entity within tolerance; every entity within tolerance at that
         * dimension (not just the nearest) is collected into the node's `geom_targets` — a node near
         * a corner shared by 2+ curves ends up with 2+ targets, by design (see `CellData.h`).
         *
         * **Edges and faces** are instead classified *topologically*, from their own boundary: an
         * edge is put on the lowest-dimensional entity containing both its corners' classifications,
         * a face on the lowest-dimensional entity containing all 4 of its edges' — see
         * `infer_targets()`. That is what keeps a cell coherent with its boundary, which proximity
         * cannot: sampling one point of an edge can otherwise land it on a curve neither of its
         * endpoints touches. Proximity remains the fallback whenever the boundary can't decide (an
         * unclassified corner, or classifications with nothing in common). **Blocks** stay pure
         * proximity, there being only volumes at that dimension.
         *
         * The single nearest target is then used to snap/refit the cell's own stored geometry: a
         * node's position is projected onto it; a curved edge's interior control points (its 2
         * endpoints stay pinned to their own, already-classified nodes) are projected onto it; a
         * face's/block's stored surface/volume is rebuilt from its (now possibly curved) boundary
         * edges/faces via `coons_surface_from_edges()`/`tfi_volume_from_faces()` — never left as a
         * stale blend of the pre-classification straight geometry.
         *
         * Cells are visited in increasing dimension because each level's inference consumes the
         * previous one's result. Not incremental: every cell is reclassified and every
         * edge/face/block geometry rebuilt from scratch on every call — see `snap_node()` for the
         * local counterpart used while editing.
         *
         * @param ATolVertex Tolerance for snapping onto a vertex — the tight threshold, since 2
         *        distinct vertices are typically much closer to each other than to any curve.
         * @param ATolCurve Tolerance for snapping onto a curve. Defaults to @p ATolVertex.
         * @param ATolSurface Tolerance for snapping onto a surface (and onto a volume). Defaults to
         *        the resolved curve tolerance.
         */
        void classify(double ATolVertex, double ATolCurve = -1.0, double ATolSurface = -1.0) {
            const Tolerances tol = resolve_tolerances(ATolVertex, ATolCurve, ATolSurface);

            for (auto it = m_cmap.template attributes<0>().begin(), itend = m_cmap.template attributes<0>().end();
                 it != itend;
                 ++it) {
                classify_node(it, tol);
            }

            for (auto it = m_cmap.template attributes<1>().begin(), itend = m_cmap.template attributes<1>().end();
                 it != itend;
                 ++it) {
                refit_edge(it, tol);
            }

            for (auto it = m_cmap.template attributes<2>().begin(), itend = m_cmap.template attributes<2>().end();
                 it != itend;
                 ++it) {
                classify_and_rebuild_face(it, tol);
            }

            for (auto it = m_cmap.template attributes<3>().begin(), itend = m_cmap.template attributes<3>().end();
                 it != itend;
                 ++it) {
                classify_and_rebuild_block(it, tol);
            }
        }

        /**
         * @brief Snaps one node onto the geometric model and brings every cell touching it back into
         * agreement — the incremental counterpart of `classify()`, meant to run when a corner is
         * released after being dragged.
         *
         * The node is reclassified by proximity and projected onto its nearest target, then every
         * incident edge, face and block is reclassified and refitted exactly as `classify()` would.
         * Because edges and faces infer their classification from their own boundary rather than
         * from a global search, re-deciding just the cells that touch @p ANode is enough: no other
         * cell's boundary changed, so no other cell's classification can have.
         *
         * @param ANode The node to snap.
         * @param ATolVertex Tolerance for snapping onto a vertex.
         * @param ATolCurve Tolerance for snapping onto a curve. Defaults to @p ATolVertex.
         * @param ATolSurface Tolerance for snapping onto a surface. Defaults to the resolved curve
         *        tolerance.
         */
        void snap_node(Node ANode, double ATolVertex, double ATolCurve = -1.0, double ATolSurface = -1.0) {
            const Tolerances tol = resolve_tolerances(ATolVertex, ATolCurve, ATolSurface);
            classify_node(ANode, tol);

            // Swept exhaustively for the same reason `move_node()` documents: a corner's own dart
            // orbit doesn't reach every cell that touches it on a still-unsewn block.
            for (auto it = m_cmap.template attributes<1>().begin(), itend = m_cmap.template attributes<1>().end();
                 it != itend;
                 ++it) {
                const Dart d = it->dart();
                const Node n0 = m_cmap.template attribute<0>(d);
                const Node n1 = m_cmap.template attribute<0>(m_cmap.template beta<1>(d));
                if (n0 == ANode || n1 == ANode) refit_edge(it, tol);
            }

            for (auto it = m_cmap.template attributes<2>().begin(), itend = m_cmap.template attributes<2>().end();
                 it != itend;
                 ++it) {
                if (face_has_node(it, ANode)) classify_and_rebuild_face(it, tol);
            }

            for (auto it = m_cmap.template attributes<3>().begin(), itend = m_cmap.template attributes<3>().end();
                 it != itend;
                 ++it) {
                if (block_has_node(it, ANode)) classify_and_rebuild_block(it, tol);
            }
        }

        /**
         * @brief Generates a full quad (2D) or hex (3D) mesh reproducing the blocking's stored
         * geometry (straight for `TEdgeCurve` degree 1, curved otherwise), subdividing every
         * edge/face/block uniformly into `ASubdivisions` intervals per parametric axis.
         *
         * Every distinct `Node`/`Edge`/`Face`/`Block` attribute is visited exactly once and its own
         * grid of sample points generated exactly once, directly from its own stored geometry
         * (`curve.value(t)` / `surface.value(u,v)` / `volume.value(u,v,w)`) — never re-derived from
         * neighbors at meshing time, and never duplicated at a seam: 2 sewn blocks share the exact
         * same `Edge`/`Face` attribute, so their common boundary's sample points are generated once
         * and reused by both sides. Which of a shared attribute's own grid indices corresponds to a
         * given block-local grid position (the "seam index") is derived purely from the corner
         * `Node` attributes' identity — exact by construction, no coordinates/tolerance involved,
         * via `bilinear_index_remap()`. A standalone 2D (quad) block becomes a mesh `Face`; a
         * genuine 3D (hex) block becomes a `Cell` — both can coexist in the returned mesh for a
         * mixed blocking. `ASubdivisions=1` reproduces exactly the coarse block topology (one
         * quad/hex per top-cell, at the block corners) — see `io::VtkMeshWriter`'s own doc comment
         * for the intended use of that as the "block structure" VTK export.
         *
         * Every mesh node also gets 2 `Variable<Int>` written under `NODE_CLASSIFICATION_DIM_VARIABLE`
         * and `NODE_CLASSIFICATION_TAG_VARIABLE`, inheriting the classification of whichever
         * Node/Edge/Face/Block attribute it was sampled from — see those constants' own doc.
         *
         * @param ASubdivisions Number of intervals to subdivide every parametric axis into; must be
         *        >= 1. Uniform across the whole blocking (a documented V1 limitation — no
         *        per-block/per-edge subdivision counts, no non-conformal adjacency).
         * @return The generated mesh.
         */
        UnstructuredMesh<CubicTraits> to_mesh(SizeT ASubdivisions) {
            assert(ASubdivisions >= 1 && "Blocking::to_mesh: ASubdivisions must be >= 1");
            const std::size_t s = ASubdivisions;
            UnstructuredMesh<CubicTraits> mesh;
            Variable<Int> &dims =
                mesh.template add_variable<Int, CellType::Node>(std::string(NODE_CLASSIFICATION_DIM_VARIABLE));
            Variable<Int> &tags =
                mesh.template add_variable<Int, CellType::Node>(std::string(NODE_CLASSIFICATION_TAG_VARIABLE));

            std::map<Node, NodeId> node_ids;
            for (auto it = m_cmap.template attributes<0>().begin(), itend = m_cmap.template attributes<0>().end();
                 it != itend;
                 ++it) {
                node_ids[it] = mesh.add_node(it->info().point);
                record_node_classification(dims, tags, node_ids[it], it->info().geom_targets);
            }

            std::map<Edge, EdgeChain> edge_chains;
            for (auto it = m_cmap.template attributes<1>().begin(), itend = m_cmap.template attributes<1>().end();
                 it != itend;
                 ++it) {
                edge_chains[it] = build_edge_chain(it, node_ids, mesh, dims, tags, s);
            }

            std::map<Face, FaceGrid> face_grids;
            for (auto it = m_cmap.template attributes<2>().begin(), itend = m_cmap.template attributes<2>().end();
                 it != itend;
                 ++it) {
                FaceGrid fg = build_face_grid(it, node_ids, edge_chains, mesh, dims, tags, s);
                if (std::find(m_hex_faces.begin(), m_hex_faces.end(), it) == m_hex_faces.end()) {
                    for (std::size_t i = 0; i < s; ++i) {
                        for (std::size_t j = 0; j < s; ++j) {
                            mesh.add_face(fg.grid[i][j], fg.grid[i + 1][j], fg.grid[i + 1][j + 1], fg.grid[i][j + 1]);
                        }
                    }
                }
                face_grids[it] = std::move(fg);
            }

            for (auto it = m_cmap.template attributes<3>().begin(), itend = m_cmap.template attributes<3>().end();
                 it != itend;
                 ++it) {
                emit_hex_cells(it, node_ids, edge_chains, face_grids, mesh, dims, tags, s);
            }

            return mesh;
        }

        /**
         * @brief Checks whether this blocking is purely 2D (contains no 3-cell/hex block) — the
         * precondition every 2D editing operation requires (this class' own `delete_face()`, and
         * the separate `FaceCollapse`/`FaceOpening`/`ChordRemoval`/`ChordInsertion` operation
         * classes built on top of `Blocking`'s public API).
         * @return true if the blocking has no 3-cells.
         */
        bool is_purely_2d() const { return nb_cells<3>() == 0; }

        /**
         * @brief Checks whether `AFace` can be deleted: currently only requires the blocking to
         * be purely 2D (see `is_purely_2d()`) — face deletion itself has no further
         * classification/topology constraint of its own (every face is unconditionally removable
         * in a dimension-3 map with no 3-cells, per CGAL's `is_removable<i>` rule for `i ==
         * dimension-1`).
         * @param AFace The face to check (unused beyond keeping the check-then-act pairing this
         *        module's other editing operations share).
         * @return true if `AFace` can be deleted.
         */
        bool can_delete_face(Face AFace) const {
            (void)AFace;
            return is_purely_2d();
        }

        /**
         * @brief Deletes `AFace` from the structure.
         *
         * Topology: a single `CGAL::Combinatorial_map::remove_cell<2>()` call. CGAL's own
         * attribute reference-counting (`Decrease_attribute_functor`, run per erased dart) already
         * garbage-collects any of `AFace`'s boundary edges/corner nodes left with zero remaining
         * incident darts as a side effect — no separate cascading-removal pass is needed here.
         * Edges/nodes still shared with another surviving face/edge are untouched and keep their
         * existing classification; a cell that *is* garbage-collected simply ceases to exist along
         * with its own classification. Geometry of every surviving cell is left exactly as it was
         * — deleting one face doesn't change what any other cell's position/curve/surface means.
         *
         * @param AFace The face to delete.
         * @pre `can_delete_face(AFace)`
         */
        void delete_face(Face AFace) {
            assert(can_delete_face(AFace) &&
                   "Blocking::delete_face: precondition violated (blocking must be purely 2D)");
            m_cmap.template remove_cell<2>(AFace->dart());
        }

        /**
         * @brief Deletes one hex block, along with every face, edge and node that existed only
         * because of it.
         *
         * Topology: a single `CGAL::Combinatorial_map::remove_cell<3>()`. Removing a cell of the
         * map's own dimension erases that cell's darts outright rather than merging it into a
         * neighbour, and CGAL's attribute reference-counting then garbage-collects whatever is left
         * holding no darts — so the block's private faces, edges and corners go with it, while
         * anything it shared with a neighbouring block survives, as a boundary cell of that
         * neighbour. Two blocks sewn along a face become one block with that face on its boundary,
         * which is what deleting one of them ought to mean.
         *
         * Unlike `delete_face()`, this has no precondition: a 3-cell in a dimension-3 map is always
         * removable, and a blocking that ends up in several disconnected pieces — or empty — is a
         * legitimate thing to be in the middle of building.
         *
         * Geometry and classification of every surviving cell are left exactly as they were:
         * deleting a block changes what exists, never what anything means.
         *
         * @param ABlock The block to delete.
         */
        void delete_block(Block ABlock) {
            // The removal renumbers whatever corner attributes it has to rebuild, which moves the
            // canonical frames derived from those ids — so the frames are noted first and the cells
            // they moved under are put back afterwards.
            run_reframed([&] { m_cmap.template remove_cell<3>(ABlock->dart()); });
            // The removal may have taken faces with it, and `m_hex_faces` would go on holding
            // handles to attributes that no longer exist — which `to_mesh()` compares against and
            // `build_connectivity()` dereferences.
            refresh_hex_faces();
        }

        // ---------------------------------------------------------------------
        // Sheet cut
        // ---------------------------------------------------------------------

        /** @brief One edge of a sheet, paired with the endpoint its cut parameter is measured from. */
        struct SheetEdge {
            /** @brief The edge itself. */
            Edge edge{};
            /** @brief The endpoint parameter 0 sits at, so one parameter cuts the whole sheet coherently. */
            Node from{};
        };

        /**
         * @brief One face a sheet crosses, with the local frame it had when the sheet was collected.
         *
         * The frame is captured up front because applying a cut invalidates it: splitting the face's
         * edges lengthens its dart cycle, so walking `beta<1>` from its dart no longer enumerates its
         * 4 corners.
         */
        struct SheetFace {
            /** @brief The face itself. */
            Face face{};
            /** @brief Its 4 corners, in its own dart-cycle order. */
            std::array<Node, 4> corners{};
            /** @brief Which of its own axes the cut runs across: 0 = local `u`, 1 = local `v`. */
            int axis = 0;
            /** @brief Local corner index of one crossing edge's `from` node, giving the cut's side. */
            int from_corner = 0;
        };

        /** @brief One block a sheet crosses, with the local frame it had when the sheet was collected.
         * @see SheetFace for why the frame is captured rather than re-walked. */
        struct SheetBlock {
            /** @brief The block itself. */
            Block block{};
            /** @brief Its 8 corners, in `HEX_CORNER_UVW` order. */
            std::array<Node, 8> corners{};
            /** @brief Which of its own axes the cut runs across: 0/1/2 = local `u`/`v`/`w`. */
            int axis = 0;
            /** @brief Local corner index of one crossing edge's `from` node, giving the cut's side. */
            int from_corner = 0;
        };

        /** @brief Everything one `cut_sheet()` would touch, as returned by `find_sheet()`. */
        struct Sheet {
            /** @brief Every edge the cut splits. */
            std::vector<SheetEdge> edges;
            /** @brief Every face the cut splits (exactly 2 of its 4 edges are in `edges`). */
            std::vector<SheetFace> faces;
            /** @brief Every block the cut splits (exactly 4 of its 12 edges are in `edges`). */
            std::vector<SheetBlock> blocks;
        };

        /**
         * @brief Collects the sheet through @p AEdge: every edge a cut of @p AEdge would have to
         * split to stay conformal, plus the faces and blocks those edges run through.
         *
         * The sheet is the transitive closure of "is parallel to, within one block" starting from
         * @p AEdge. Since 2 sewn blocks share the very same `Edge` attribute, that relation walks
         * from block to block on its own — which is what makes the cut propagate through the whole
         * layer rather than stopping at the selected block. A standalone quad block joins in the
         * same way, through its 2 parallel edges instead of a hex's 4.
         *
         * Each edge is collected together with the endpoint its parameter must be measured from, so
         * that one parameter cuts every edge on the same side. That side is carried across blocks by
         * integer corner coordinates (`HEX_CORNER_UVW`/`QUAD_CORNER_IJ`), never by geometry, so a
         * neighbour block whose local frame runs the other way still gets `1-t` rather than `t`
         * without any distance or tolerance being involved.
         *
         * @param AEdge The edge the cut is aimed at; parameter 0 is its own stored curve's start
         *        (`curve.control_points()[0]`), matching what `cut_sheet()`'s parameter means.
         * @return The sheet, or `std::nullopt` when it is not homogeneously cuttable — either a face
         *         carries a number of sheet edges other than 0 or 2, a block other than 0 or 4, or
         *         one edge is reached from 2 blocks that disagree on which end it starts at. All 3
         *         mean the sheet closes back onto itself in a way that leaves no single well-defined
         *         cut, and are reported rather than cut wrongly.
         */
        std::optional<Sheet> find_sheet(Edge AEdge) {
            std::map<Edge, Node> start_of;
            std::vector<Edge> pending;
            start_of.emplace(AEdge, curve_start_node(AEdge));
            pending.push_back(AEdge);

            while (!pending.empty()) {
                const Edge e = pending.back();
                pending.pop_back();
                const Node from = start_of.at(e);

                for (const Block b : blocks_of_edge(e)) {
                    const std::array<Node, 8> corners = block_local_nodes(b);
                    const std::array<Edge, 12> edges = block_edges(b, corners);
                    const std::size_t ie = table_index(edges, e);
                    const int axis = hex_edge_axis(ie);
                    const int side = HEX_CORNER_UVW[static_cast<std::size_t>(node_index(corners, from))][axis];
                    for (std::size_t k = 0; k < 12; ++k) {
                        if (hex_edge_axis(k) != axis) continue;
                        const auto a = static_cast<std::size_t>(HEX_EDGES[k].first);
                        const auto b2 = static_cast<std::size_t>(HEX_EDGES[k].second);
                        const Node ns = (HEX_CORNER_UVW[a][axis] == side) ? corners[a] : corners[b2];
                        if (!extend_sheet(start_of, pending, edges[k], ns)) return std::nullopt;
                    }
                }

                for (const Face f : standalone_faces_of_edge(e)) {
                    const std::array<Node, 4> corners = face_local_nodes(f);
                    const std::array<Edge, 4> edges = face_edges(f, corners);
                    const std::size_t ie = table_index(edges, e);
                    const int axis = quad_edge_axis(ie);
                    const int side = QUAD_CORNER_IJ[static_cast<std::size_t>(node_index(corners, from))][axis];
                    for (std::size_t k = 0; k < 4; ++k) {
                        if (quad_edge_axis(k) != axis) continue;
                        const auto a = static_cast<std::size_t>(QUAD_EDGES[k].first);
                        const auto b2 = static_cast<std::size_t>(QUAD_EDGES[k].second);
                        const Node ns = (QUAD_CORNER_IJ[a][axis] == side) ? corners[a] : corners[b2];
                        if (!extend_sheet(start_of, pending, edges[k], ns)) return std::nullopt;
                    }
                }
            }

            Sheet sheet;
            sheet.edges.reserve(start_of.size());
            for (const auto &[e, from] : start_of) {
                sheet.edges.push_back(SheetEdge{e, from});
            }
            // Ordered by the ids of the nodes each edge joins, never by the order `start_of` happens
            // to hold them in: that map is keyed on attribute handles, handles are addresses, and
            // addresses move between runs. The order matters because the cut hands out node ids in
            // it, and the canonical frames are built on those ids — so leaving it to the allocator
            // makes the whole operation answer differently from one run to the next.
            std::sort(sheet.edges.begin(), sheet.edges.end(), [this](const SheetEdge &a, const SheetEdge &b) {
                return edge_order(a.edge) < edge_order(b.edge);
            });
            if (!collect_sheet_faces(start_of, sheet) || !collect_sheet_blocks(start_of, sheet)) {
                return std::nullopt;
            }
            return sheet;
        }

        /**
         * @brief Measures the volume one block encloses, from its own stored geometry.
         *
         * The block's `volume.value(u,v,w)` is sampled on an `(S+1)^3` grid and the resulting `S^3`
         * sub-hexes are each split into 6 tetrahedra around a main diagonal, whose signed volumes are
         * summed. That makes the answer **exact** for a block whose faces are planar — an axis-aligned
         * box at `ASubdivisions = 1` measures its extent exactly — and an approximation that converges
         * as `ASubdivisions` grows for a curved or warped one, whose bilinear faces no tetrahedral
         * split can follow exactly.
         *
         * The result is signed, so a block whose frame ended up inverted reports a negative volume
         * rather than silently reading as valid.
         *
         * @param ABlock The block to measure.
         * @param ASubdivisions Intervals per parametric axis; must be >= 1.
         * @return Its signed volume.
         */
        double block_volume(Block ABlock, SizeT ASubdivisions = 1) {
            assert(ASubdivisions >= 1 && "Blocking::block_volume: ASubdivisions must be >= 1");
            const std::size_t s = ASubdivisions;
            const auto at = [&](std::size_t i, std::size_t j, std::size_t k) {
                const double d = static_cast<double>(s);
                return ABlock->info().volume.value(
                    static_cast<double>(i) / d, static_cast<double>(j) / d, static_cast<double>(k) / d);
            };

            // The 6 tetrahedra around the 0-6 diagonal, which tile a hexahedron of planar faces.
            static constexpr std::array<std::array<std::size_t, 4>, 6> TETS = {std::array<std::size_t, 4>{0, 1, 2, 6},
                                                                               {0, 2, 3, 6},
                                                                               {0, 3, 7, 6},
                                                                               {0, 7, 4, 6},
                                                                               {0, 4, 5, 6},
                                                                               {0, 5, 1, 6}};
            double total = 0.0;
            for (std::size_t i = 0; i < s; ++i) {
                for (std::size_t j = 0; j < s; ++j) {
                    for (std::size_t k = 0; k < s; ++k) {
                        // The sub-hex's 8 corners, in the same HEX8 order HEX_CORNER_UVW describes.
                        const std::array<Point3d, 8> c = {at(i, j, k),
                                                          at(i + 1, j, k),
                                                          at(i + 1, j + 1, k),
                                                          at(i, j + 1, k),
                                                          at(i, j, k + 1),
                                                          at(i + 1, j, k + 1),
                                                          at(i + 1, j + 1, k + 1),
                                                          at(i, j + 1, k + 1)};
                        for (const auto &t : TETS) {
                            const Vector3d u(c[t[0]], c[t[1]]);
                            const Vector3d v(c[t[0]], c[t[2]]);
                            const Vector3d w(c[t[0]], c[t[3]]);
                            total += u.cross(v).dot(w) / 6.0;
                        }
                    }
                }
            }
            return total;
        }

        /**
         * @brief Measures every block of the blocking, in the order `attributes<3>()` traverses them.
         * @param ASubdivisions Intervals per parametric axis; must be >= 1 (see `block_volume()`).
         * @return One signed volume per block.
         */
        std::vector<double> block_volumes(SizeT ASubdivisions = 1) {
            std::vector<double> volumes;
            for (auto it = m_cmap.template attributes<3>().begin(), itend = m_cmap.template attributes<3>().end();
                 it != itend;
                 ++it) {
                volumes.push_back(block_volume(it, ASubdivisions));
            }
            return volumes;
        }

        /**
         * @brief Where a cut at @p AParam falls along one sheet edge, as a parameter along that
         * edge's *own* stored curve.
         *
         * A sheet measures its parameter from `SheetEdge::from`, but each edge's stored curve runs
         * whichever way it was built — so on an edge whose curve starts at the far end, the same cut
         * sits at `1 - AParam`. Every part of the cut goes through here, which is what keeps the
         * whole sheet cut on one side.
         *
         * @param AEdge One edge of a sheet, from `find_sheet()`.
         * @param AParam The cut parameter, as `cut_sheet()` defines it.
         * @return The matching parameter along @p AEdge's own curve.
         */
        double cut_parameter(const SheetEdge &AEdge, double AParam) {
            return (curve_start_node(AEdge.edge) == AEdge.from) ? AParam : 1.0 - AParam;
        }

        /**
         * @brief The point a cut at @p AParam would fall on, along one sheet edge.
         * @param AEdge One edge of a sheet, from `find_sheet()`.
         * @param AParam The cut parameter, as `cut_sheet()` defines it.
         * @return That point, on the edge's own curve.
         */
        Point3d cut_point(const SheetEdge &AEdge, double AParam) {
            return AEdge.edge->info().curve.value(cut_parameter(AEdge, AParam));
        }

        /**
         * @brief Cuts the blocking along the whole sheet through @p AEdge, at parameter @p AParam.
         *
         * Every edge of the sheet gains a node, every face it crosses an edge, every block it crosses
         * a face — so one hex becomes 2, and the cut propagates to every block sharing those edges
         * (see `find_sheet()`). The result stays conformal: 2 blocks that were sewn are still sewn,
         * through the 2 halves of the face they shared.
         *
         * **Curvature is preserved exactly.** Each half keeps the parent's own geometry, restricted:
         * an edge's curve, a face's surface and a block's volume are De Casteljau *subdivisions* of
         * what they were cut out of (see `BezierCurve::split()`), never a refit from the new
         * boundaries. Sampling the cut blocking therefore traces the very same points as sampling it
         * before — the point at `(u,v,w)` of a block cut at `s` along `u` is the point at `(u/s,v,w)`
         * of its low half. Re-deriving the halves through `coons_surface_from_edges()` instead would
         * move the surface, since the restriction of a Coons patch to a sub-rectangle is not the
         * Coons patch of that sub-rectangle's boundary.
         *
         * Classification is inherited rather than recomputed, which needs no geometric query and
         * cannot invent a target: each half of a cell keeps the cell's own `geom_targets`, and a cell
         * born inside another (the node on an edge, the edge inside a face, the face inside a block)
         * takes its parent's.
         *
         * Worth knowing: a later `classify()` rebuilds every face's surface and block's volume from
         * their boundaries, which is that method's documented job — and which discards the exact
         * subdivision geometry this operation was careful to keep. Cut first, classify after, is the
         * order that keeps both.
         *
         * @param AEdge The edge the cut is aimed at.
         * @param AParam Where along @p AEdge to cut, measured along its own stored curve, strictly
         *        inside `(0, 1)`. Every other edge of the sheet is cut at the matching parameter from
         *        the matching side, so the cut stays homogeneous across blocks whose local frames
         *        disagree.
         * @return false, changing nothing at all, if @p AParam is not strictly inside `(0, 1)` or the
         *         sheet is not homogeneously cuttable (see `find_sheet()`).
         */
        bool cut_sheet(Edge AEdge, double AParam) {
            if (!(AParam > 0.0 && AParam < 1.0)) return false;
            const std::optional<Sheet> sheet = find_sheet(AEdge);
            if (!sheet.has_value()) return false;

            // Ordered by increasing dimension, and each pass finished before the next starts: a
            // face's own frame stops being walkable the moment one of its edges is split, so every
            // frame the cut needs was captured by find_sheet() before any of this ran.
            // A cut renumbers corner attributes too — CGAL rebuilds whatever orbit its insertions
            // disturb — so the frames of the cells it is *not* reshaping are noted and put back, the
            // same way a deletion does it.
            // Each phase gets its own record-and-restore. One pass around the whole cut is not
            // enough: a phase legitimately reshapes some cells and stores their geometry in the frame
            // they have at that moment, and a *later* phase can then renumber a corner under them.
            // Wrapped phase by phase, every cell is either the one being reshaped — skipped, its
            // geometry having just been written — or one merely renumbered, and put back.
            std::map<Node, MidNode> mids;
            run_reframed([&] { split_sheet_edges(*sheet, AParam, mids); });
            run_reframed([&] { split_sheet_faces(*sheet, AParam, mids); });
            run_reframed([&] { split_sheet_blocks(*sheet, AParam, mids); });
            // A cut both splits faces in 2 and creates one per block, so which faces belong to a
            // block has changed; re-derived rather than tracked through each of those steps.
            refresh_hex_faces();
            return true;
        }

        /**
         * @brief Moves one node to a new position, refitting the geometry of every cell it belongs
         * to so the structure stays consistent.
         *
         * Every incident edge's curve is re-derived as a straight interpolation between its 2
         * endpoints' current positions (interior control points included, for a curved
         * `TEdgeCurve`), then every incident face's surface and block's volume is rebuilt from
         * those edges. Classification (`geom_targets`) is deliberately left untouched throughout:
         * moving a node is a purely geometric edit, and re-deciding what a cell is classified on is
         * `classify()`'s job — call it afterwards to re-snap onto the geometric model.
         *
         * @param ANode The node to move.
         * @param ANewPosition Its new position.
         */
        void move_node(Node ANode, const Point3d &ANewPosition) {
            ANode->info().point = ANewPosition;

            // Swept exhaustively rather than through `one_dart_per_incident_cell<i,0>`: a vertex's
            // own dart orbit only reaches its other incident cells through beta2/beta3, so on a
            // standalone (still unsewn) block a corner's orbit is a single dart and the incidence
            // iterators silently miss most of the cells that actually touch it. Same full-sweep
            // shape as `classify()`; a blocking is small enough for that to be irrelevant.
            for (auto it = m_cmap.template attributes<1>().begin(), itend = m_cmap.template attributes<1>().end();
                 it != itend;
                 ++it) {
                const Dart d = it->dart();
                const Node n0 = m_cmap.template attribute<0>(d);
                const Node n1 = m_cmap.template attribute<0>(m_cmap.template beta<1>(d));
                if (n0 != ANode && n1 != ANode) continue;
                store_curve(it, straight_curve(n0->info().point, n1->info().point), n0);
            }

            for (auto it = m_cmap.template attributes<2>().begin(), itend = m_cmap.template attributes<2>().end();
                 it != itend;
                 ++it) {
                if (face_has_node(it, ANode)) rebuild_face_surface(it);
            }

            for (auto it = m_cmap.template attributes<3>().begin(), itend = m_cmap.template attributes<3>().end();
                 it != itend;
                 ++it) {
                if (block_has_node(it, ANode)) rebuild_block_volume(it);
            }
        }

        /** @brief Gives access to the underlying combinatorial map. @return The internal map. */
        Map &cmap() { return m_cmap; }
        /** @copydoc cmap() */
        const Map &cmap() const { return m_cmap; }

        /** @brief Gives access to the associated geometric model. @return The geometric model. */
        const TGeomModel &geom_model() const { return *m_geom_model; }

        /**
         * @brief Gets the number of `TDim`-cells in the blocking structure.
         * @tparam TDim The dimension of cells to count, in [0,3].
         * @return The number of `TDim`-cells.
         */
        template<unsigned int TDim>
        std::size_t nb_cells() const {
            return m_cmap.template number_of_attributes<TDim>();
        }

        /** @brief Checks the topological validity of the underlying map. @return true if valid. */
        bool is_valid_topology() const { return m_cmap.is_valid(); }

    private:
        /** @brief The (u,v) grid coordinates (unscaled, in `{0,1}`) of a quad's 4 local corners
         * 0..3, matching `QUAD_EDGES`' `u` along corner 0->1 and `v` along corner 0->3. */
        static constexpr std::array<std::array<int, 2>, 4> QUAD_CORNER_IJ = {std::array{0, 0},
                                                                             std::array{1, 0},
                                                                             std::array{1, 1},
                                                                             std::array{0, 1}};

        /** @brief The (u,v,w) grid coordinates (unscaled, in `{0,1}`) of a hex's 8 local corners
         * 0..7 — derived from `CubicTraits::Cell`'s HEX8 layout and the existing, tested
         * `volume.value(u,v,w)` corner correspondences (see `blocking_creation_tests.cpp`). */
        static constexpr std::array<std::array<int, 3>, 8> HEX_CORNER_UVW = {std::array{0, 0, 0},
                                                                             std::array{1, 0, 0},
                                                                             std::array{1, 1, 0},
                                                                             std::array{0, 1, 0},
                                                                             std::array{0, 0, 1},
                                                                             std::array{1, 0, 1},
                                                                             std::array{1, 1, 1},
                                                                             std::array{0, 1, 1}};

        /** @brief A `to_mesh()`-generated node-id chain along one edge attribute's own curve,
         * parameter `t=i/S` at index `i`, plus which physical `Node` corresponds to index 0 — since
         * an edge's own curve direction isn't necessarily aligned with any one particular
         * face's/block's local traversal of it (see `oriented_curve()`).
         */
        struct EdgeChain {
            /** @brief The node at parameter t=0 (`ids[0]`). */
            Node start{};
            /** @brief One mesh node id per sample point, index `i` at parameter `t=i/S`. */
            std::vector<NodeId> ids;
        };

        /** @brief A `to_mesh()`-generated `(S+1)x(S+1)` grid of node ids for one face attribute,
         * indexed in the face's own local `(u,v)` frame (see `classify_and_rebuild_face()`), plus
         * the 4 local corner nodes that frame is anchored to (needed to reconcile this frame against
         * a differently-rotated/reflected view of the same face from an incident block).
         */
        struct FaceGrid {
            /** @brief The face's own 4 local corner nodes, walked from its own dart. */
            std::array<Node, 4> local_nodes{};
            /** @brief `(S+1)x(S+1)` node ids, `grid[i][j]` at parameter `(u,v)=(i/S,j/S)`. */
            std::vector<std::vector<NodeId>> grid;
        };

        /** @brief A `to_mesh()`-generated `(S+1)^3` grid of node ids for one hex block. */
        using Grid3D = std::vector<std::vector<std::vector<NodeId>>>;

        /**
         * @brief Sews every pair of still-free, `AK`-corner-matching `TDim`-cells among @p ACells
         * (position-based matching: 2 blocks built independently at the same location have
         * *separate*, coincident-position attributes before sewing).
         * @tparam TDim The sew dimension (3 for hex faces, matching `AK=4` corners each).
         * @param ACells The candidate cell attributes (all of the same, already-known dimension).
         * @param AK Number of corners each cell has (4 for a hex's bounding face).
         */
        template<unsigned int TDim>
        void sew_matching(std::vector<Face> &ACells, int AK) {
            for (std::size_t i = 0; i < ACells.size(); ++i) {
                const Dart dA = ACells[i]->dart();
                if (!m_cmap.template is_free<TDim>(dA)) continue;
                for (std::size_t j = i + 1; j < ACells.size(); ++j) {
                    if (!m_cmap.template is_free<TDim>(ACells[j]->dart())) continue;
                    Dart dB = ACells[j]->dart();
                    bool matched = false;
                    for (int rot = 0; rot < AK; ++rot) {
                        if (positions_match_reversed(dA, dB, AK)) {
                            matched = true;
                            break;
                        }
                        dB = m_cmap.template beta<1>(dB);
                    }
                    if (matched) {
                        m_cmap.template sew<TDim>(dA, dB);
                        break;
                    }
                }
            }
        }

        /**
         * @brief Checks whether 2 free (unsewn) edges represent the same segment traversed in
         * opposite directions — the correct 2-sewability condition: `dA`'s start matches `dB`'s
         * end, and `dA`'s end matches `dB`'s start.
         *
         * Deliberately *not* `positions_match_reversed()`: that helper (correct for `sew_matching`'s
         * `AK=4` hex-face case, which searches all `AK` rotations of `dB` around its own face until
         * one aligns) compares `dA`'s start against `dB`'s *own* start, then `dA`'s end against the
         * start of the dart *before* `dB` in `dB`'s face — a check whose 2 comparisons only land on
         * the true pair of segment endpoints for a specific rotation of `dB`. A free edge has just
         * one dart (no "other rotation" to try the way a face's own boundary does — walking
         * `beta<1>` off of it moves to a *different* edge in the same face entirely, not another
         * orientation of this one), so there's no rotation search that fixes it: the correct check
         * has to compare directly against `dB`'s actual end/start. Using `positions_match_reversed`
         * here was verified to silently mis-sew 2 quads sharing an edge — gluing one quad's
         * unrelated edge to the other's — while `is_valid_topology()` still reported the corrupted
         * map as valid (topological validity alone doesn't catch a wrong *geometric*
         * correspondence).
         */
        bool free_edges_match_reversed(Dart ADartA, Dart ADartB) {
            return m_cmap.template attribute<0>(ADartA)->info().point ==
                       m_cmap.template attribute<0>(m_cmap.template beta<1>(ADartB))->info().point &&
                   m_cmap.template attribute<0>(m_cmap.template beta<1>(ADartA))->info().point ==
                       m_cmap.template attribute<0>(ADartB)->info().point;
        }

        /**
         * @brief Sews every pair of still-free, position-matching edges among all standalone quad
         * blocks' boundary edges. A hex's own 12 edges are always already 2-sewn (to their 2
         * sibling faces within that hex, done at creation time), so `is_free<2>` naturally selects
         * only standalone quad edges here — no separate bookkeeping needed to exclude hex edges.
         */
        void sew_free_quad_edges() {
            std::vector<Dart> free_edges;
            for (auto it = m_cmap.template attributes<1>().begin(), itend = m_cmap.template attributes<1>().end();
                 it != itend;
                 ++it) {
                Dart d = it->dart();
                if (m_cmap.template is_free<2>(d)) free_edges.push_back(d);
            }
            for (std::size_t i = 0; i < free_edges.size(); ++i) {
                const Dart dA = free_edges[i];
                if (!m_cmap.template is_free<2>(dA)) continue;
                for (std::size_t j = i + 1; j < free_edges.size(); ++j) {
                    const Dart dB = free_edges[j];
                    if (!m_cmap.template is_free<2>(dB)) continue;
                    if (free_edges_match_reversed(dA, dB)) {
                        m_cmap.template sew<2>(dA, dB);
                        break;
                    }
                }
            }
        }

        /**
         * @brief Checks whether walking `AK` corners forward from @p ADartA (via `beta<1>`) matches
         * walking `AK` corners forward from @p ADartB but in reverse (via `beta<0>`) — the
         * orientation 2 abutting cells' shared boundary must have for a geometrically-correct sew
         * (they're traversed in opposite cyclic order from either side).
         * @param ADartA First cell's starting dart.
         * @param ADartB Second cell's starting dart.
         * @param AK Number of corners to check.
         * @return true if all `AK` corner positions match in this reversed correspondence.
         */
        bool positions_match_reversed(Dart ADartA, Dart ADartB, int AK) {
            Dart a = ADartA;
            Dart b = ADartB;
            for (int s = 0; s < AK; ++s) {
                if (!(m_cmap.template attribute<0>(a)->info().point == m_cmap.template attribute<0>(b)->info().point)) {
                    return false;
                }
                a = m_cmap.template beta<1>(a);
                b = m_cmap.template beta<0>(b);
            }
            return true;
        }

        /**
         * @brief Resolves `classify()`/`snap_node()`'s defaulted tolerance arguments: an omitted
         * (negative) curve tolerance falls back to the vertex one, and an omitted surface tolerance
         * to the curve one.
         * @param ATolVertex Vertex tolerance.
         * @param ATolCurve Curve tolerance, or negative to reuse @p ATolVertex.
         * @param ATolSurface Surface tolerance, or negative to reuse the resolved curve tolerance.
         * @return The resolved tolerances.
         */
        static Tolerances resolve_tolerances(double ATolVertex, double ATolCurve, double ATolSurface) {
            const double curve = (ATolCurve < 0.0) ? ATolVertex : ATolCurve;
            const double surface = (ATolSurface < 0.0) ? curve : ATolSurface;
            return Tolerances{ATolVertex, curve, surface};
        }

        /**
         * @brief Outcome of searching `m_geom_model` for entities near a query point: every entity
         * found within tolerance (see `classify_position()`), plus the dimension/tag of the single
         * nearest one (used to actually snap/project geometry).
         */
        struct ClassifyResult {
            /** @brief Every entity found within tolerance, as (dimension, entity_tag) pairs. */
            std::vector<std::pair<GroupDim, Int>> targets;
            /** @brief Dimension of the nearest entity in `targets`, or `Undefined` if `targets` is empty. */
            GroupDim nearest_dim = GroupDim::Undefined;
            /** @brief `entity_tag()` of the nearest entity in `targets`. */
            Int nearest_tag{};
            /** @brief Checks whether any entity was found. @return true if `targets` is non-empty. */
            [[nodiscard]] bool any() const { return nearest_dim != GroupDim::Undefined; }
        };

        /**
         * @brief Scans one dimension's worth of `m_geom_model` entities, collecting every one within
         * @p ATol of @p AP into @p AResult and updating its nearest-entity bookkeeping.
         * @tparam TEntityRange A range of `GeomEntityConcept`-conforming entities.
         * @param AEntities The entities to scan (e.g. `m_geom_model->curves()`).
         * @param ADim The dimension @p AEntities belongs to.
         * @param AP The query point.
         * @param ATol Distance tolerance.
         * @param AResult Accumulates matching entities and the running nearest one.
         * @param ABestDist Running best (smallest) distance seen so far across calls for this result.
         */
        template<typename TEntityRange>
        static void scan_dimension(const TEntityRange &AEntities,
                                   GroupDim ADim,
                                   const Point3d &AP,
                                   double ATol,
                                   ClassifyResult &AResult,
                                   double &ABestDist) {
            for (const auto &entity : AEntities) {
                const double d = entity.distance(AP);
                if (d <= ATol) {
                    AResult.targets.emplace_back(ADim, entity.entity_tag());
                    if (!AResult.any() || d < ABestDist) {
                        ABestDist = d;
                        AResult.nearest_dim = ADim;
                        AResult.nearest_tag = entity.entity_tag();
                    }
                }
            }
        }

        /**
         * @brief Searches `m_geom_model` for entities of dimension `>= AMinDim` near @p AP, dimension
         * by dimension, stopping at the lowest dimension with any match — see `classify()`.
         * @param AMinDim Lowest entity dimension to consider (the classified cell's own topological
         *        dimension).
         * @param AP The query point.
         * @param ATol The per-dimension snapping tolerances.
         * @return The classification outcome (possibly empty, if nothing is within tolerance anywhere).
         */
        ClassifyResult classify_position(GroupDim AMinDim, const Point3d &AP, const Tolerances &ATol) const {
            ClassifyResult result;
            double best = 0.0;
            if (AMinDim <= GroupDim::Dim0) {
                scan_dimension(m_geom_model->vertices(), GroupDim::Dim0, AP, ATol.vertex, result, best);
                if (result.any()) return result;
            }
            if (AMinDim <= GroupDim::Dim1) {
                scan_dimension(m_geom_model->curves(), GroupDim::Dim1, AP, ATol.curve, result, best);
                if (result.any()) return result;
            }
            if (AMinDim <= GroupDim::Dim2) {
                scan_dimension(m_geom_model->surfaces(), GroupDim::Dim2, AP, ATol.surface, result, best);
                if (result.any()) return result;
            }
            if (AMinDim <= GroupDim::Dim3) {
                // Volumes reuse the surface tolerance: a `GeomModelConcept` volume's `distance()` is
                // free to report 0 everywhere (as `FacetedVolume`'s stub does), which makes any
                // threshold inert here — so there is nothing a 4th tolerance could usefully control.
                scan_dimension(m_geom_model->volumes(), GroupDim::Dim3, AP, ATol.surface, result, best);
            }
            return result;
        }

        /**
         * @brief Classifies one node by proximity and projects it onto its nearest target, leaving
         * it where it is when nothing is within tolerance.
         * @param ANode The node to classify.
         * @param ATol The per-dimension snapping tolerances.
         */
        void classify_node(Node ANode, const Tolerances &ATol) {
            const auto result = classify_position(GroupDim::Dim0, ANode->info().point, ATol);
            ANode->info().geom_targets = result.targets;
            if (result.any()) {
                project_onto(result.nearest_dim, result.nearest_tag, ANode->info().point);
            }
        }

        /**
         * @brief Gathers every geometric entity containing *any* of @p ATargets — the union of the
         * model's `containing_entities()` over a cell's classification, which is the set of entities
         * that cell could still be said to lie on.
         * @param ATargets One cell's `geom_targets`.
         * @return The union, as a sorted set.
         */
        std::set<std::pair<GroupDim, Int>> containing_set(const std::vector<std::pair<GroupDim, Int>> &ATargets) const {
            std::set<std::pair<GroupDim, Int>> result;
            for (const auto &[dim, tag] : ATargets) {
                for (const auto &entry : m_geom_model->containing_entities(dim, tag)) {
                    result.insert(entry);
                }
            }
            return result;
        }

        /**
         * @brief Infers what a cell is classified on from the classification of its own boundary:
         * the lowest-dimensional geometric entities containing *every* boundary cell's own
         * classification.
         *
         * This is what keeps a cell's classification coherent with its boundary's, which proximity
         * alone cannot guarantee — sampling a single point of an edge can put it on a curve its own
         * two endpoints have nothing to do with. An edge whose 2 corners sit on the same geometric
         * curve belongs on that curve; one whose corners sit on 2 different curves bounding a common
         * surface belongs on that surface. Intersecting the boundary's containing sets (see
         * `containing_set()`) states exactly that, and the model computes those sets from shared
         * backing mesh cells, so no tolerance enters here at all.
         *
         * Returns *every* candidate at the winning dimension rather than one, matching how
         * `classify_position()` reports several equally-plausible targets; picking the single
         * nearest of them for projection is `nearest_of()`'s job.
         *
         * @param ABoundary Each boundary cell's `geom_targets`. An empty entry (an unclassified
         *        boundary cell) makes inference impossible, since it constrains nothing.
         * @param AMinDim The cell's own topological dimension — it can never be classified below it.
         * @return The inferred targets, or empty when inference is impossible (an unclassified
         *         boundary cell, or boundary classifications with nothing in common), in which case
         *         the caller falls back to a proximity search.
         */
        std::vector<std::pair<GroupDim, Int>>
        infer_targets(const std::vector<const std::vector<std::pair<GroupDim, Int>> *> &ABoundary,
                      GroupDim AMinDim) const {
            if (ABoundary.empty()) return {};

            std::set<std::pair<GroupDim, Int>> common;
            bool first = true;
            for (const auto *targets : ABoundary) {
                if (targets->empty()) return {}; // Unclassified boundary: nothing to infer from.
                const auto candidates = containing_set(*targets);
                if (first) {
                    common = candidates;
                    first = false;
                } else {
                    std::set<std::pair<GroupDim, Int>> intersection;
                    std::ranges::set_intersection(common, candidates, std::inserter(intersection, intersection.end()));
                    common = std::move(intersection);
                }
                if (common.empty()) return {};
            }

            // `common` is sorted by (dimension, tag), so the first entry at or above the cell's own
            // dimension already sits at the winning — lowest — dimension.
            GroupDim best = GroupDim::Undefined;
            std::vector<std::pair<GroupDim, Int>> result;
            for (const auto &[dim, tag] : common) {
                if (dim < AMinDim) continue;
                if (best == GroupDim::Undefined) best = dim;
                if (dim != best) break;
                result.emplace_back(dim, tag);
            }
            return result;
        }

        /**
         * @brief Picks the entity of @p ATargets closest to @p AP — the geometric tie-break applied
         * when `infer_targets()` leaves several equally-valid candidates at the winning dimension
         * (e.g. an edge lying on a curve shared by 2 surfaces), and the entity a cell's geometry is
         * then projected onto.
         * @param ATargets Candidate entities.
         * @param AP Representative point of the cell.
         * @return The nearest candidate, or `Undefined` if @p ATargets is empty.
         */
        ClassifyResult nearest_of(const std::vector<std::pair<GroupDim, Int>> &ATargets, const Point3d &AP) const {
            ClassifyResult result;
            result.targets = ATargets;
            double best = 0.0;
            for (const auto &[dim, tag] : ATargets) {
                const double d = distance_to(dim, tag, AP);
                if (!result.any() || d < best) {
                    best = d;
                    result.nearest_dim = dim;
                    result.nearest_tag = tag;
                }
            }
            return result;
        }

        /**
         * @brief Distance from @p AP to the geometric entity identified by (@p ADim, @p ATag).
         * @param ADim Dimension of the entity.
         * @param ATag `entity_tag()` of the entity.
         * @param AP The query point.
         * @return The distance, or infinity if no such entity exists.
         */
        double distance_to(GroupDim ADim, Int ATag, const Point3d &AP) const {
            if (ADim == GroupDim::Dim0) {
                if (const auto *v = m_geom_model->vertex_by_tag(ATag)) return v->distance(AP);
            } else if (ADim == GroupDim::Dim1) {
                if (const auto *c = m_geom_model->curve_by_tag(ATag)) return c->distance(AP);
            } else if (ADim == GroupDim::Dim2) {
                if (const auto *s = m_geom_model->surface_by_tag(ATag)) return s->distance(AP);
            } else if (ADim == GroupDim::Dim3) {
                if (const auto *vol = m_geom_model->volume_by_tag(ATag)) return vol->distance(AP);
            }
            return std::numeric_limits<double>::infinity();
        }

        /**
         * @brief Projects @p AP onto the geometric entity identified by (@p ADim, @p ATag), in place.
         * @param ADim Dimension of the entity to project onto (see `GeomModelConcept`'s `*_by_tag()`).
         * @param ATag `entity_tag()` of the entity to project onto.
         * @param AP The point to snap, modified in place. Left unchanged if no such entity exists.
         */
        void project_onto(GroupDim ADim, Int ATag, Point3d &AP) const {
            if (ADim == GroupDim::Dim0) {
                if (const auto *v = m_geom_model->vertex_by_tag(ATag)) v->project(AP);
            } else if (ADim == GroupDim::Dim1) {
                if (const auto *c = m_geom_model->curve_by_tag(ATag)) c->project(AP);
            } else if (ADim == GroupDim::Dim2) {
                if (const auto *s = m_geom_model->surface_by_tag(ATag)) s->project(AP);
            } else if (ADim == GroupDim::Dim3) {
                if (const auto *vol = m_geom_model->volume_by_tag(ATag)) vol->project(AP);
            }
        }

        /**
         * @brief Classifies one edge and refits its curve so that the curve itself lies on what it
         * is classified on, its 2 endpoints pinned to their (already-classified) nodes.
         *
         * Classification is inferred from the edge's own 2 corner nodes via `infer_targets()`,
         * falling back to a proximity search at its midpoint when they can't decide (an
         * unclassified corner, or corners with no common containing entity).
         *
         * The fit **interpolates points sampled on the geometry** rather than projecting the control
         * points onto it. Projecting the control points is the tempting shortcut, and it is wrong: a
         * Bezier curve passes through its 2 endpoints but *not* through its interior control points,
         * staying strictly inside their convex hull — so control points sitting exactly on a curved
         * geometric curve leave the actual edge bowed off it, always falling short of the geometry
         * it is supposed to follow. Instead, points on the geometry are obtained first (by projecting
         * a straight-line guess), and the control points that make the curve pass *through* those
         * points are then solved for — see `interpolating_curve()`.
         *
         * @param AEdge The edge to classify and refit.
         * @param ATol The per-dimension snapping tolerances, used only by the fallback search.
         */
        void refit_edge(Edge AEdge, const Tolerances &ATol) {
            const Dart d = AEdge->dart();
            const Node n0 = m_cmap.template attribute<0>(d);
            const Node n1 = m_cmap.template attribute<0>(m_cmap.template beta<1>(d));
            const Point3d &p0 = n0->info().point;
            const Point3d &p1 = n1->info().point;

            const Point3d midpoint = AEdge->info().curve.value(0.5);
            const auto inferred = infer_targets({&n0->info().geom_targets, &n1->info().geom_targets}, GroupDim::Dim1);
            const auto result =
                inferred.empty() ? classify_position(GroupDim::Dim1, midpoint, ATol) : nearest_of(inferred, midpoint);
            AEdge->info().geom_targets = result.targets;

            const std::size_t n = m_degree + 1;
            const Vector3d chord(p0, p1);

            // The points the fitted curve must pass through: its 2 endpoints, plus interior samples
            // pulled onto the geometry. Without a classification there is nothing to pull them onto,
            // so they stay on the chord and the fit collapses to the straight edge.
            std::vector<Point3d> samples(n);
            for (std::size_t i = 0; i < n; ++i) {
                const double t = (n > 1) ? static_cast<double>(i) / static_cast<double>(n - 1) : 0.0;
                Point3d s = p0 + chord * t;
                if (result.any() && i > 0 && i + 1 < n) {
                    project_onto(result.nearest_dim, result.nearest_tag, s);
                }
                samples[i] = s;
            }
            samples[0] = p0;
            samples[n - 1] = p1;

            // On a curve the geometry also dictates which *direction* the edge must leave its ends
            // in, and honouring that matters as much as passing through the right points — a fit
            // free to choose its end tangents picks badly wrong ones (~30 degrees out on a plain
            // circular arc). On a surface there is no single such direction, so plain interpolation
            // stands.
            if (result.nearest_dim == GroupDim::Dim1) {
                const auto *curve = m_geom_model->curve_by_tag(result.nearest_tag);
                if (curve != nullptr) {
                    // Far more samples than there are unknowns, and parameterized by their own
                    // spacing rather than uniformly: a point projected from the chord at parameter
                    // `t` does not sit at parameter `t` along the geometry, and assuming it does
                    // biases the fit (it leaves the tangents pointing the right way but too short).
                    std::vector<Point3d> dense;
                    dense.reserve(4 * n);
                    for (std::size_t i = 0; i < 4 * n; ++i) {
                        const double t = static_cast<double>(i) / static_cast<double>(4 * n - 1);
                        Point3d s = p0 + chord * t;
                        project_onto(result.nearest_dim, result.nearest_tag, s);
                        dense.push_back(s);
                    }
                    dense.front() = p0;
                    dense.back() = p1;

                    store_curve(AEdge,
                                tangent_constrained_curve(m_degree,
                                                          dense,
                                                          chord_length_parameters(dense),
                                                          oriented_tangent(curve->tangent(p0), chord),
                                                          oriented_tangent(curve->tangent(p1), chord)),
                                n0);
                    return;
                }
            }
            store_curve(AEdge, interpolating_curve(samples), n0);
        }

        /**
         * @brief Assigns each of @p APoints a curve parameter in `[0,1]` proportional to the
         * distance travelled along them — the standard chord-length parameterization, which matches
         * where points actually fall along a curve far better than spacing them uniformly.
         * @param APoints The sample points, in order.
         * @return One parameter per point, starting at 0 and ending at 1. All zero for a degenerate
         *         (zero-length) sample run.
         */
        static std::vector<double> chord_length_parameters(const std::vector<Point3d> &APoints) {
            std::vector<double> parameters(APoints.size(), 0.0);
            for (std::size_t i = 1; i < APoints.size(); ++i) {
                parameters[i] = parameters[i - 1] + Vector3d(APoints[i - 1], APoints[i]).norm();
            }
            const double total = parameters.back();
            if (total <= 0.0) return parameters;
            for (double &p : parameters) {
                p /= total;
            }
            return parameters;
        }

        /**
         * @brief Flips @p ATangent, whose sign is arbitrary (see `FacetedCurve::tangent()`), to
         * travel the same way as @p AChord.
         * @param ATangent The tangent to orient.
         * @param AChord The direction of travel along the edge.
         * @return The oriented tangent.
         */
        static Vector3d oriented_tangent(const Vector3d &ATangent, const Vector3d &AChord) {
            return (ATangent.dot(AChord) < 0.0) ? -ATangent : ATangent;
        }

        /**
         * @brief Builds the curve leaving @p ASamples' 2 endpoints along @p AStartTangent and
         * @p AEndTangent, fitted as closely as possible to the interior samples.
         *
         * The end tangents are imposed exactly by construction: a Bezier's start tangent is
         * `P1 - P0` and its end tangent `Pn - P(n-1)`, so pinning `P1 = P0 + a * AStartTangent` and
         * `P(n-1) = Pn - b * AEndTangent` fixes both directions whatever `a` and `b` turn out to be.
         * Those 2 lengths are then the only freedom left (for degree 3; higher degrees keep their
         * extra middle control points free too), and are chosen by least squares against the
         * interior samples — the classic tangent-constrained Bezier fit.
         *
         * @param ADegree The degree to fit at.
         * @param ASamples The points to fit, first and last being the endpoints.
         * @param AParameters The curve parameter each sample should be reproduced at, one per entry
         *        of @p ASamples — see `chord_length_parameters()`.
         * @param AStartTangent Unit direction the curve must leave `ASamples.front()` in.
         * @param AEndTangent Unit direction the curve must arrive at `ASamples.back()` in.
         * @return The fitted curve. Falls back to plain interpolation if either tangent is null (a
         *         degenerate faceting) or the fit is singular.
         */
        static TEdgeCurve tangent_constrained_curve(std::size_t ADegree,
                                                    const std::vector<Point3d> &ASamples,
                                                    const std::vector<double> &AParameters,
                                                    const Vector3d &AStartTangent,
                                                    const Vector3d &AEndTangent) {
            const std::size_t n = ADegree + 1;
            const std::size_t degree = ADegree;
            std::vector<Point3d> ends(n);
            ends[0] = ASamples.front();
            ends[n - 1] = ASamples.back();
            if (n < 3) {
                return interpolating_curve(ends);
            } else {
                if (AStartTangent.norm_sq() < 1e-24 || AEndTangent.norm_sq() < 1e-24) {
                    return interpolating_curve(ends);
                }

                const Point3d &p0 = ASamples.front();
                const Point3d &pn = ASamples.back();

                // Least squares over the interior samples for the 2 tangent lengths. Each sample
                // contributes its 3 coordinates, so even degree 3 (2 unknowns) is overdetermined.
                double ata[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
                double atb[2] = {0.0, 0.0};
                for (std::size_t i = 1; i + 1 < ASamples.size(); ++i) {
                    const double t = AParameters[i];
                    // B(t) = base(t) + a * B_1(t) * start - b * B_(n-2)(t) * end, with everything
                    // that does not depend on a or b gathered into base(t).
                    const double ca = bernstein(degree, 1, t);
                    const double cb = -bernstein(degree, degree - 1, t);
                    Vector3d base = Vector3d(Point3d(0, 0, 0), p0) * (bernstein(degree, 0, t) + ca) +
                                    Vector3d(Point3d(0, 0, 0), pn) * (bernstein(degree, degree, t) - cb);
                    // Degree 4 and up keep middle control points of their own; leaving them on the
                    // chord keeps this a 2-unknown fit while the ends stay exact.
                    for (std::size_t j = 2; j + 2 < n; ++j) {
                        const double s = static_cast<double>(j) / static_cast<double>(degree);
                        const Point3d mid = p0 + Vector3d(p0, pn) * s;
                        base += Vector3d(Point3d(0, 0, 0), mid) * bernstein(degree, j, t);
                    }

                    const Vector3d residual = Vector3d(Point3d(0, 0, 0), ASamples[i]) - base;
                    const Vector3d da = AStartTangent * ca;
                    const Vector3d db = AEndTangent * cb;
                    ata[0][0] += da.dot(da);
                    ata[0][1] += da.dot(db);
                    ata[1][0] += db.dot(da);
                    ata[1][1] += db.dot(db);
                    atb[0] += da.dot(residual);
                    atb[1] += db.dot(residual);
                }

                const double det = ata[0][0] * ata[1][1] - ata[0][1] * ata[1][0];
                if (std::abs(det) < 1e-18) {
                    return interpolating_curve(ends);
                }
                const double a = (atb[0] * ata[1][1] - ata[0][1] * atb[1]) / det;
                const double b = (ata[0][0] * atb[1] - atb[0] * ata[1][0]) / det;

                TEdgeCurve fitted(degree);
                fitted[0] = p0;
                fitted[n - 1] = pn;
                fitted[1] = p0 + AStartTangent * a;
                fitted[n - 2] = pn - AEndTangent * b;
                for (std::size_t j = 2; j + 2 < n; ++j) {
                    const double s = static_cast<double>(j) / static_cast<double>(degree);
                    fitted[j] = p0 + Vector3d(p0, pn) * s;
                }
                return fitted;
            }
        }

        /**
         * @brief Builds the curve of degree `TEdgeCurve::Degree` passing through @p ASamples at
         * uniform parameters `t_i = i/(n-1)`.
         *
         * Solves `sum_j B_j(t_i) * P_j = S_i` for the control points `P`, where `B_j` are the
         * Bernstein basis polynomials. The 2 endpoints are already control points (a Bezier
         * interpolates them), so only the interior ones are unknown, leaving an `(n-2)x(n-2)` system
         * — at most 3x3 here, small enough for plain Gaussian elimination with partial pivoting and
         * not worth pulling in a linear algebra dependency for.
         *
         * @param ASamples The `n` points to pass through, first and last being the endpoints.
         * @return The interpolating curve. Returns the samples as-is when there is no interior
         *         control point to solve for (degree 1), where the curve already interpolates them.
         */
        static TEdgeCurve interpolating_curve(const std::vector<Point3d> &ASamples) {
            const std::size_t n = ASamples.size();
            const std::size_t degree = n - 1;
            const std::size_t unknowns = (n > 2) ? n - 2 : 0;

            TEdgeCurve fitted(degree);
            fitted[0] = ASamples[0];
            fitted[n - 1] = ASamples[n - 1];
            if (unknowns == 0) {
                return fitted;
            } else {
                // Augmented system: `unknowns` equations, one per interior sample, with the 3
                // coordinates carried as 3 right-hand sides (the matrix is the same for all of them).
                std::vector<std::vector<double>> sys(unknowns, std::vector<double>(unknowns + 3, 0.0));
                for (std::size_t r = 0; r < unknowns; ++r) {
                    const double t = static_cast<double>(r + 1) / static_cast<double>(degree);
                    for (std::size_t c = 0; c < unknowns; ++c) {
                        sys[r][c] = bernstein(degree, c + 1, t);
                    }
                    // Move the 2 known endpoint terms over to the right-hand side.
                    const Vector3d rhs = Vector3d(Point3d(0, 0, 0), ASamples[r + 1]) -
                                         Vector3d(Point3d(0, 0, 0), ASamples[0]) * bernstein(degree, 0, t) -
                                         Vector3d(Point3d(0, 0, 0), ASamples[n - 1]) * bernstein(degree, degree, t);
                    sys[r][unknowns] = rhs.x();
                    sys[r][unknowns + 1] = rhs.y();
                    sys[r][unknowns + 2] = rhs.z();
                }

                for (std::size_t col = 0; col < unknowns; ++col) {
                    std::size_t pivot = col;
                    for (std::size_t r = col + 1; r < unknowns; ++r) {
                        if (std::abs(sys[r][col]) > std::abs(sys[pivot][col])) pivot = r;
                    }
                    std::swap(sys[col], sys[pivot]);
                    // A Bernstein collocation matrix at distinct parameters is non-singular, so this
                    // only guards against a degenerate build rather than an expected case.
                    if (std::abs(sys[col][col]) < 1e-300) continue;
                    for (std::size_t r = 0; r < unknowns; ++r) {
                        if (r == col) continue;
                        const double factor = sys[r][col] / sys[col][col];
                        for (std::size_t c = col; c < unknowns + 3; ++c) {
                            sys[r][c] -= factor * sys[col][c];
                        }
                    }
                }

                for (std::size_t r = 0; r < unknowns; ++r) {
                    const double d = sys[r][r];
                    fitted[r + 1] = Point3d(sys[r][unknowns] / d, sys[r][unknowns + 1] / d, sys[r][unknowns + 2] / d);
                }
                return fitted;
            }
        }

        /**
         * @brief Evaluates the Bernstein basis polynomial `B_i^n(t) = C(n,i) t^i (1-t)^(n-i)`.
         * @param ADegree The basis degree `n`.
         * @param AIndex The basis index `i`, in `[0, ADegree]`.
         * @param AT The parameter.
         * @return The value of the polynomial.
         */
        static double bernstein(std::size_t ADegree, std::size_t AIndex, double AT) {
            double binomial = 1.0;
            for (std::size_t k = 0; k < AIndex; ++k) {
                binomial = binomial * static_cast<double>(ADegree - k) / static_cast<double>(k + 1);
            }
            return binomial * std::pow(AT, static_cast<double>(AIndex)) *
                   std::pow(1.0 - AT, static_cast<double>(ADegree - AIndex));
        }

        /**
         * @brief Builds a control-point-reversed copy of a curve.
         * @param ACurve The curve to reverse.
         * @return The reversed curve.
         */
        static TEdgeCurve reversed_curve(const TEdgeCurve &ACurve) {
            TEdgeCurve rev(ACurve.degree());
            const std::size_t n = ACurve.nb_control_points();
            for (std::size_t i = 0; i < n; ++i) {
                rev[i] = ACurve[n - 1 - i];
            }
            return rev;
        }

        /**
         * @brief Returns @p AEdge's curve oriented to start at @p AExpectedStart, reversing it when
         * the edge stores it the other way round — used to orient an already-fitted edge curve
         * consistently for `coons_surface_from_edges()`/`tfi_volume_from_faces()` regardless of which
         * of its 2 incident faces/blocks is asking.
         *
         * The comparison is on node *identity*, never on position. Which way an edge stores its curve
         * is a combinatorial fact — `store_curve()` keeps it pinned to the edge's own dart — and
         * asking geometry about it instead is both slower and wrong: 2 ways of computing the same
         * corner agree as real numbers and not as doubles, so an equality test on coordinates
         * silently answers "the other end" for any cell whose curve and nodes were built by different
         * routes (see `curve_start_node()`).
         *
         * @param AEdge The edge whose curve is wanted.
         * @param AExpectedStart The node its first control point must sit on.
         * @return That curve, oriented to start at @p AExpectedStart.
         */
        TEdgeCurve oriented_curve(Edge AEdge, Node AExpectedStart) {
            if (curve_start_node(AEdge) == AExpectedStart) return AEdge->info().curve;
            return reversed_curve(AEdge->info().curve);
        }

        /**
         * @brief The node an edge's stored curve starts at.
         *
         * Purely combinatorial, and deliberately anchored to the *lower of the edge's 2 endpoint
         * handles* rather than to its own dart: an edge carries darts pointing both ways, which of
         * them the attribute names is CGAL's business, and CGAL does re-point it — during the very
         * cell insertions a cut is made of. Ordering the 2 endpoints instead depends on nothing that
         * moves. Attribute handles already order this way elsewhere here (they key the maps in
         * `to_mesh()` and `find_sheet()`).
         *
         * Position is never consulted. Asking geometry which end a curve starts at is what broke
         * repeated cuts: an edge born inside a face takes its curve from a surface's control grid
         * while its nodes were placed by evaluating another curve, so the 2 agree as real numbers and
         * not as doubles, and an equality test on coordinates then answers "the other end".
         *
         * @param AEdge The edge to inspect.
         * @return Its curve's start node.
         */
        Node curve_start_node(Edge AEdge) {
            const Dart d = AEdge->dart();
            const Node a = m_cmap.template attribute<0>(d);
            const Node b = m_cmap.template attribute<0>(m_cmap.template beta<1>(d));
            // Which dart is asked decides only *which pair* comes back, never its order — every dart
            // of an edge spans the same 2 nodes.
            return (a < b) ? a : b;
        }

        /**
         * @brief Stores a curve on an edge, oriented to the class' invariant — reversing @p ACurve
         * when @p AStart is not the end `curve_start_node()` names.
         *
         * Every write of an edge's curve goes through here, which is what lets `curve_start_node()`
         * answer by looking rather than by measuring.
         *
         * @param AEdge The edge to write to.
         * @param ACurve The curve, running from @p AStart to the edge's other endpoint.
         * @param AStart The node @p ACurve starts at.
         */
        void store_curve(Edge AEdge, const TEdgeCurve &ACurve, Node AStart) {
            AEdge->info().curve = (curve_start_node(AEdge) == AStart) ? ACurve : reversed_curve(ACurve);
        }

        /**
         * @brief Classifies one face and rebuilds its stored surface via `coons_surface_from_edges()`
         * from its 4 (now possibly refitted) boundary edges. Applies uniformly to a standalone 2D
         * (quad) block's own face and to one of a 3D (hex) block's 6 bounding faces — both are
         * ordinary 2-cells with exactly 4 boundary edges.
         *
         * Classification is inferred from those 4 edges' own classifications via `infer_targets()`
         * (so a face all of whose edges lie on one surface lands on that surface), falling back to a
         * proximity search at the face's corner centroid when they can't decide. Edges must
         * therefore already be classified — `classify()` and `snap_node()` both do faces after
         * edges for that reason.
         *
         * @param AFace The face to classify and rebuild.
         * @param ATol The per-dimension snapping tolerances, used only by the fallback search.
         */
        void classify_and_rebuild_face(Face AFace, const Tolerances &ATol) {
            const Dart fd = AFace->dart();
            std::array<Node, 4> local_nodes{};
            Dart walk = fd;
            for (std::size_t c = 0; c < 4; ++c) {
                local_nodes[c] = m_cmap.template attribute<0>(walk);
                walk = m_cmap.template beta<1>(walk);
            }

            Vector3d acc(0.0, 0.0, 0.0);
            for (const Node &node : local_nodes) {
                acc += Vector3d(local_nodes[0]->info().point, node->info().point);
            }
            const Point3d center = local_nodes[0]->info().point + acc * 0.25;

            std::vector<const std::vector<std::pair<GroupDim, Int>> *> boundary;
            for (auto it = m_cmap.template one_dart_per_incident_cell<1, 2>(fd).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<1, 2>(fd).end();
                 it != itend;
                 ++it) {
                boundary.push_back(&m_cmap.template attribute<1>(it)->info().geom_targets);
            }

            const auto inferred = infer_targets(boundary, GroupDim::Dim2);
            const auto result =
                inferred.empty() ? classify_position(GroupDim::Dim2, center, ATol) : nearest_of(inferred, center);
            AFace->info().geom_targets = result.targets;

            rebuild_face_surface(AFace);
        }

        /**
         * @brief Checks whether @p ANode is one of @p AFace's 4 corners.
         * @param AFace The face to inspect.
         * @param ANode The node to look for.
         * @return true if the face has that corner.
         */
        bool face_has_node(Face AFace, Node ANode) {
            Dart walk = AFace->dart();
            for (std::size_t c = 0; c < 4; ++c) {
                if (m_cmap.template attribute<0>(walk) == ANode) return true;
                walk = m_cmap.template beta<1>(walk);
            }
            return false;
        }

        /**
         * @brief Checks whether @p ANode is one of @p ABlock's 8 corners.
         * @param ABlock The block to inspect.
         * @param ANode The node to look for.
         * @return true if the block has that corner.
         */
        bool block_has_node(Block ABlock, Node ANode) {
            const Dart bd = ABlock->dart();
            for (auto it = m_cmap.template one_dart_per_incident_cell<0, 3>(bd).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<0, 3>(bd).end();
                 it != itend;
                 ++it) {
                if (m_cmap.template attribute<0>(it) == ANode) return true;
            }
            return false;
        }

        /**
         * @brief Rebuilds one face's stored surface via `coons_surface_from_edges()` from its 4
         * boundary edges' current curves, leaving its classification untouched — the pure-geometry
         * half of `classify_and_rebuild_face()`, also used on its own by `move_node()`.
         * @param AFace The face to rebuild.
         */
        void rebuild_face_surface(Face AFace) {
            const Dart fd = AFace->dart();
            // The face's own frame, canonical rather than whichever dart happens to name it — see
            // `face_local_nodes()`. Every reader and writer of a face's surface goes through it, or
            // one of them ends up reading in a frame another one did not write.
            const std::array<Node, 4> local_nodes = face_local_nodes(AFace);

            std::array<TEdgeCurve, 4> curves{};
            for (auto it = m_cmap.template one_dart_per_incident_cell<1, 2>(fd).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<1, 2>(fd).end();
                 it != itend;
                 ++it) {
                const int a = node_index(local_nodes, m_cmap.template attribute<0>(it));
                const int b = node_index(local_nodes, m_cmap.template attribute<0>(m_cmap.template beta<1>(it)));
                const std::size_t k = find_edge(QUAD_EDGES, a, b);
                const Edge e = m_cmap.template attribute<1>(it);
                curves[k] = oriented_curve(e, local_nodes[static_cast<std::size_t>(QUAD_EDGES[k].first)]);
            }
            AFace->info().surface = coons_surface_from_edges(curves[0], curves[2], curves[3], curves[1]);
        }

        /**
         * @brief Classifies one block against `m_geom_model` (dimension >= 3 only, using the block's
         * own corner centroid) and rebuilds its stored volume via `tfi_volume_from_faces()`. The 6
         * bounding face surfaces it needs are re-derived directly from the block's own (now possibly
         * refitted) 12 boundary edges, in the block's own local frame — deliberately not read back
         * from the boundary `Face` attributes' own stored surfaces, which may have been rebuilt by
         * `classify_and_rebuild_face()` in an independent local orientation (reconciling the two is
         * `to_mesh()`'s seam-indexing concern, not this method's).
         * Classification stays a plain proximity search here, unlike edges' and faces': the only
         * entities at dimension 3 are volumes, and inferring from the block's 6 faces could only
         * ever return those same volumes.
         *
         * @param ABlock The block to classify and rebuild.
         * @param ATol The per-dimension snapping tolerances.
         */
        void classify_and_rebuild_block(Block ABlock, const Tolerances &ATol) {
            const Dart bd = ABlock->dart();
            std::array<Node, 8> local_nodes{};
            local_nodes[0] = m_cmap.template attribute<0>(bd);
            local_nodes[1] = m_cmap.template attribute<0>(m_cmap.template beta<1>(bd));
            local_nodes[2] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1>(bd));
            local_nodes[3] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1, 1>(bd));
            local_nodes[4] = m_cmap.template attribute<0>(m_cmap.template beta<2, 1, 1>(bd));
            local_nodes[5] = m_cmap.template attribute<0>(m_cmap.template beta<1, 2, 1, 1>(bd));
            local_nodes[6] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1, 2, 1, 1>(bd));
            local_nodes[7] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1, 1, 2, 1, 1>(bd));

            Vector3d acc(0.0, 0.0, 0.0);
            for (const Node &node : local_nodes) {
                acc += Vector3d(local_nodes[0]->info().point, node->info().point);
            }
            const Point3d center = local_nodes[0]->info().point + acc * (1.0 / 8.0);
            const auto result = classify_position(GroupDim::Dim3, center, ATol);
            ABlock->info().geom_targets = result.targets;

            rebuild_block_volume(ABlock);
        }

        /**
         * @brief Rebuilds one block's stored volume via `tfi_volume_from_faces()` from its 12
         * boundary edges' current curves, leaving its classification untouched — the pure-geometry
         * half of `classify_and_rebuild_block()`, also used on its own by `move_node()`.
         * @param ABlock The block to rebuild.
         */
        void rebuild_block_volume(Block ABlock) {
            const Dart bd = ABlock->dart();
            std::array<Node, 8> local_nodes{};
            local_nodes[0] = m_cmap.template attribute<0>(bd);
            local_nodes[1] = m_cmap.template attribute<0>(m_cmap.template beta<1>(bd));
            local_nodes[2] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1>(bd));
            local_nodes[3] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1, 1>(bd));
            local_nodes[4] = m_cmap.template attribute<0>(m_cmap.template beta<2, 1, 1>(bd));
            local_nodes[5] = m_cmap.template attribute<0>(m_cmap.template beta<1, 2, 1, 1>(bd));
            local_nodes[6] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1, 2, 1, 1>(bd));
            local_nodes[7] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1, 1, 2, 1, 1>(bd));

            std::array<TEdgeCurve, 12> curves{};
            for (auto it = m_cmap.template one_dart_per_incident_cell<1, 3>(bd).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<1, 3>(bd).end();
                 it != itend;
                 ++it) {
                const int a = node_index(local_nodes, m_cmap.template attribute<0>(it));
                const int b = node_index(local_nodes, m_cmap.template attribute<0>(m_cmap.template beta<1>(it)));
                const std::size_t k = find_edge(HEX_EDGES, a, b);
                const Edge e = m_cmap.template attribute<1>(it);
                curves[k] = oriented_curve(e, local_nodes[static_cast<std::size_t>(HEX_EDGES[k].first)]);
            }

            std::array<FaceSurfaceT, 6> face_surfaces{};
            for (std::size_t f = 0; f < 6; ++f) {
                const auto &spec = HEX_FACES[f];
                face_surfaces[f] = coons_surface_from_edges(
                    curves[spec.edge_u0], curves[spec.edge_u1], curves[spec.edge_v0], curves[spec.edge_v1]);
            }
            ABlock->info().volume = tfi_volume_from_faces(face_surfaces[0],
                                                          face_surfaces[1],
                                                          face_surfaces[2],
                                                          face_surfaces[3],
                                                          face_surfaces[4],
                                                          face_surfaces[5]);
        }

        /**
         * @brief Writes one mesh node's classification into `ADims`/`ATags`: the dimension and tag
         * of @p ATargets' first entry (see `NODE_CLASSIFICATION_DIM_VARIABLE`'s doc for the "first
         * target speaks for all" convention, also used by `BlockingFacade`'s `classification_dim()`),
         * or -1/-1 when @p ATargets is empty (unclassified).
         * @param ADims The `NODE_CLASSIFICATION_DIM_VARIABLE` variable, modified in place.
         * @param ATags The `NODE_CLASSIFICATION_TAG_VARIABLE` variable, modified in place.
         * @param AId The mesh node to record classification for.
         * @param ATargets The owning block-structure cell's `geom_targets`.
         */
        static void record_node_classification(Variable<Int> &ADims,
                                               Variable<Int> &ATags,
                                               NodeId AId,
                                               const std::vector<std::pair<GroupDim, Int>> &ATargets) {
            if (ATargets.empty()) {
                ADims[AId.value] = -1;
                ATags[AId.value] = -1;
            } else {
                ADims[AId.value] = static_cast<Int>(ATargets.front().first);
                ATags[AId.value] = ATargets.front().second;
            }
        }

        /**
         * @brief Builds one edge's `to_mesh()` node-id chain: its 2 endpoints reuse the already
         * mapped corner node ids, its `S-1` interior points are newly added mesh nodes sampled
         * directly from the edge's own curve, classified on `AEdge` itself.
         * @param AEdge The edge to sample.
         * @param ANodeIds Already-populated node-attribute-to-mesh-id map (its corner endpoints).
         * @param AMesh The mesh being built, appended to for interior points.
         * @param ADims The `NODE_CLASSIFICATION_DIM_VARIABLE` variable, extended for interior points.
         * @param ATags The `NODE_CLASSIFICATION_TAG_VARIABLE` variable, extended for interior points.
         * @param AS Subdivisions per axis.
         * @return The built chain.
         */
        EdgeChain build_edge_chain(Edge AEdge,
                                   const std::map<Node, NodeId> &ANodeIds,
                                   UnstructuredMesh<CubicTraits> &AMesh,
                                   Variable<Int> &ADims,
                                   Variable<Int> &ATags,
                                   std::size_t AS) {
            const Dart d = AEdge->dart();
            const Node na = m_cmap.template attribute<0>(d);
            const Node nb = m_cmap.template attribute<0>(m_cmap.template beta<1>(d));
            // Asked of the class' own orientation invariant rather than by comparing the curve's
            // first control point against a node's coordinates: the 2 can differ in their last bits
            // for any edge whose curve and nodes were built by different routes, and the chain would
            // then be laid down backwards.
            const Node start = curve_start_node(AEdge);

            EdgeChain chain;
            chain.start = start;
            const Node end = (start == na) ? nb : na;
            chain.ids.resize(AS + 1);
            chain.ids[0] = ANodeIds.at(chain.start);
            chain.ids[AS] = ANodeIds.at(end);
            for (std::size_t i = 1; i < AS; ++i) {
                const double t = static_cast<double>(i) / static_cast<double>(AS);
                chain.ids[i] = AMesh.add_node(AEdge->info().curve.value(t));
                record_node_classification(ADims, ATags, chain.ids[i], AEdge->info().geom_targets);
            }
            return chain;
        }

        /**
         * @brief Builds one face's `to_mesh()` `(S+1)x(S+1)` node-id grid: its 4 corners and 4
         * boundary rows/columns reuse already-built node ids (corners, and the appropriate
         * already-built `EdgeChain`, oriented to match this face's own local frame); only its
         * `(S-1)^2` strictly interior points are newly added mesh nodes, sampled from the face's own
         * surface.
         * @param AFace The face to sample.
         * @param ANodeIds Already-populated node-attribute-to-mesh-id map.
         * @param AEdgeChains Already-populated edge-attribute-to-chain map (its 4 boundary edges).
         * @param AMesh The mesh being built, appended to for interior points.
         * @param ADims The `NODE_CLASSIFICATION_DIM_VARIABLE` variable, extended for interior points.
         * @param ATags The `NODE_CLASSIFICATION_TAG_VARIABLE` variable, extended for interior points.
         * @param AS Subdivisions per axis.
         * @return The built grid.
         */
        FaceGrid build_face_grid(Face AFace,
                                 const std::map<Node, NodeId> &ANodeIds,
                                 const std::map<Edge, EdgeChain> &AEdgeChains,
                                 UnstructuredMesh<CubicTraits> &AMesh,
                                 Variable<Int> &ADims,
                                 Variable<Int> &ATags,
                                 std::size_t AS) {
            FaceGrid fg;
            const Dart fd = AFace->dart();
            // Same canonical frame the surface was stored in — see `face_local_nodes()`.
            fg.local_nodes = face_local_nodes(AFace);

            fg.grid.assign(AS + 1, std::vector<NodeId>(AS + 1));
            fg.grid[0][0] = ANodeIds.at(fg.local_nodes[0]);
            fg.grid[AS][0] = ANodeIds.at(fg.local_nodes[1]);
            fg.grid[AS][AS] = ANodeIds.at(fg.local_nodes[2]);
            fg.grid[0][AS] = ANodeIds.at(fg.local_nodes[3]);

            for (auto it = m_cmap.template one_dart_per_incident_cell<1, 2>(fd).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<1, 2>(fd).end();
                 it != itend;
                 ++it) {
                const int a = node_index(fg.local_nodes, m_cmap.template attribute<0>(it));
                const int b = node_index(fg.local_nodes, m_cmap.template attribute<0>(m_cmap.template beta<1>(it)));
                const std::size_t k = find_edge(QUAD_EDGES, a, b);
                const Edge e = m_cmap.template attribute<1>(it);
                fill_boundary_line(fg.grid, QUAD_EDGES[k], fg.local_nodes, AEdgeChains.at(e), AS);
            }

            for (std::size_t i = 1; i < AS; ++i) {
                for (std::size_t j = 1; j < AS; ++j) {
                    const double u = static_cast<double>(i) / static_cast<double>(AS);
                    const double v = static_cast<double>(j) / static_cast<double>(AS);
                    fg.grid[i][j] = AMesh.add_node(AFace->info().surface.value(u, v));
                    record_node_classification(ADims, ATags, fg.grid[i][j], AFace->info().geom_targets);
                }
            }
            return fg;
        }

        /**
         * @brief Fills one boundary row/column of a face's `(S+1)x(S+1)` grid from an edge's chain,
         * oriented to match this face's local `(u,v)` frame.
         * @param AGrid The face's grid, modified in place.
         * @param AEdge The `QUAD_EDGES` entry (local corner index pair) for this boundary edge.
         * @param ALocalNodes The face's own 4 local corner nodes.
         * @param AChain The edge's own node-id chain.
         * @param AS Subdivisions per axis.
         */
        static void fill_boundary_line(std::vector<std::vector<NodeId>> &AGrid,
                                       const std::pair<int, int> &AEdge,
                                       const std::array<Node, 4> &ALocalNodes,
                                       const EdgeChain &AChain,
                                       std::size_t AS) {
            const auto a = static_cast<std::size_t>(AEdge.first);
            const auto b = static_cast<std::size_t>(AEdge.second);
            const int varying_axis = (QUAD_CORNER_IJ[a][0] != QUAD_CORNER_IJ[b][0]) ? 0 : 1;
            const auto fixed_axis = static_cast<std::size_t>(1 - varying_axis);
            const auto fixed_value = static_cast<std::size_t>(QUAD_CORNER_IJ[a][fixed_axis]);
            const Node zero_end =
                (QUAD_CORNER_IJ[a][static_cast<std::size_t>(varying_axis)] == 0) ? ALocalNodes[a] : ALocalNodes[b];
            for (std::size_t m = 0; m <= AS; ++m) {
                const NodeId value = (AChain.start == zero_end) ? AChain.ids[m] : AChain.ids[AS - m];
                if (varying_axis == 0) {
                    AGrid[m][fixed_value * AS] = value;
                } else {
                    AGrid[fixed_value * AS][m] = value;
                }
            }
        }

        /**
         * @brief Fills one edge-line of a hex block's `(S+1)^3` grid from an edge's chain, oriented
         * to match the block's own local `(u,v,w)` frame.
         * @param AGrid The block's grid, modified in place.
         * @param AEdge The `HEX_EDGES` entry (local corner index pair) for this boundary edge.
         * @param ALocalNodes The block's own 8 local corner nodes.
         * @param AChain The edge's own node-id chain.
         * @param AS Subdivisions per axis.
         */
        static void fill_hex_edge_line(Grid3D &AGrid,
                                       const std::pair<int, int> &AEdge,
                                       const std::array<Node, 8> &ALocalNodes,
                                       const EdgeChain &AChain,
                                       std::size_t AS) {
            const auto a = static_cast<std::size_t>(AEdge.first);
            const auto b = static_cast<std::size_t>(AEdge.second);
            std::size_t varying_axis = 0;
            for (std::size_t ax = 0; ax < 3; ++ax) {
                if (HEX_CORNER_UVW[a][ax] != HEX_CORNER_UVW[b][ax]) {
                    varying_axis = ax;
                    break;
                }
            }
            std::array<std::size_t, 2> fixed_axes{};
            std::size_t fc = 0;
            for (std::size_t ax = 0; ax < 3; ++ax) {
                if (ax != varying_axis) fixed_axes[fc++] = ax;
            }
            const auto fv0 = static_cast<std::size_t>(HEX_CORNER_UVW[a][fixed_axes[0]]);
            const auto fv1 = static_cast<std::size_t>(HEX_CORNER_UVW[a][fixed_axes[1]]);
            const Node zero_end = (HEX_CORNER_UVW[a][varying_axis] == 0) ? ALocalNodes[a] : ALocalNodes[b];

            for (std::size_t m = 0; m <= AS; ++m) {
                const NodeId value = (AChain.start == zero_end) ? AChain.ids[m] : AChain.ids[AS - m];
                std::array<std::size_t, 3> idx{};
                idx[varying_axis] = m;
                idx[fixed_axes[0]] = fv0 * AS;
                idx[fixed_axes[1]] = fv1 * AS;
                AGrid[idx[0]][idx[1]][idx[2]] = value;
            }
        }

        /**
         * @brief Remaps a `(i,j)` grid position from one square's local `(0..S)^2` index frame to
         * another's, given the 4 target-frame coordinates (scaled to `0`/`AS`) that the source
         * frame's own 4 local corners 0..3 (at `(0,0)`,`(AS,0)`,`(AS,AS)`,`(0,AS)`) map to — an exact
         * affine (bilinear-on-a-square) remap, since corner coordinates are always exactly `0`/`AS`.
         * Used both to place a hex block's own `(u,v,w)` grid index for a boundary face's interior
         * point, and to look up that point's node id in the face's own independently-oriented grid —
         * see `fill_hex_face_interior()`.
         * @param ACornerTarget The target-frame `(p,q)` for source corners 0,1,2,3 respectively.
         * @param Ai Source frame's first (`u`-like) index, in `[0,AS]`.
         * @param Aj Source frame's second (`v`-like) index, in `[0,AS]`.
         * @param AS Subdivisions per axis.
         * @return The remapped `(p,q)` in the target frame.
         */
        static std::pair<std::size_t, std::size_t>
        bilinear_index_remap(const std::array<std::pair<long, long>, 4> &ACornerTarget,
                             std::size_t Ai,
                             std::size_t Aj,
                             std::size_t AS) {
            const auto s_l = static_cast<long>(AS);
            const auto i_l = static_cast<long>(Ai);
            const auto j_l = static_cast<long>(Aj);
            const auto [p0, q0] = ACornerTarget[0];
            const auto [p1, q1] = ACornerTarget[1];
            const auto [p3, q3] = ACornerTarget[3];
            const long p = p0 + ((p1 - p0) / s_l) * i_l + ((p3 - p0) / s_l) * j_l;
            const long q = q0 + ((q1 - q0) / s_l) * i_l + ((q3 - q0) / s_l) * j_l;
            return {static_cast<std::size_t>(p), static_cast<std::size_t>(q)};
        }

        /**
         * @brief Fills one boundary face's strictly-interior `(S-1)^2` points of a hex block's
         * `(S+1)^3` grid, reusing the face's own already-built `FaceGrid` (no new mesh nodes) via
         * `bilinear_index_remap()` for both the block-local `(u,v,w)` placement and the face-local
         * lookup.
         * @param AGrid The block's grid, modified in place.
         * @param ACornerIdx The face's 4 local corner indices (into `ALocalNodes`), in the order
         *        walked from the block's own dart for this face.
         * @param ALocalNodes The block's own 8 local corner nodes.
         * @param AFaceGrid The face's own already-built grid.
         * @param AS Subdivisions per axis.
         */
        static void fill_hex_face_interior(Grid3D &AGrid,
                                           const std::array<int, 4> &ACornerIdx,
                                           const std::array<Node, 8> &ALocalNodes,
                                           const FaceGrid &AFaceGrid,
                                           std::size_t AS) {
            std::size_t fixed_axis = 0;
            for (std::size_t ax = 0; ax < 3; ++ax) {
                const int v0 = HEX_CORNER_UVW[static_cast<std::size_t>(ACornerIdx[0])][ax];
                bool all_same = true;
                for (int c : ACornerIdx) {
                    if (HEX_CORNER_UVW[static_cast<std::size_t>(c)][ax] != v0) {
                        all_same = false;
                        break;
                    }
                }
                if (all_same) {
                    fixed_axis = ax;
                    break;
                }
            }
            const auto fixed_value =
                static_cast<std::size_t>(HEX_CORNER_UVW[static_cast<std::size_t>(ACornerIdx[0])][fixed_axis]);
            std::array<std::size_t, 2> varying_axes{};
            std::size_t vc = 0;
            for (std::size_t ax = 0; ax < 3; ++ax) {
                if (ax != fixed_axis) varying_axes[vc++] = ax;
            }

            static constexpr std::array<std::pair<long, long>, 4> FACE_CORNER_PQ = {std::pair<long, long>{0, 0},
                                                                                    std::pair<long, long>{1, 0},
                                                                                    std::pair<long, long>{1, 1},
                                                                                    std::pair<long, long>{0, 1}};

            std::array<std::pair<long, long>, 4> block_to_uvw{};
            std::array<std::pair<long, long>, 4> block_to_face_pq{};
            for (std::size_t k = 0; k < 4; ++k) {
                const auto c = static_cast<std::size_t>(ACornerIdx[k]);
                block_to_uvw[k] = {static_cast<long>(HEX_CORNER_UVW[c][varying_axes[0]]) * static_cast<long>(AS),
                                   static_cast<long>(HEX_CORNER_UVW[c][varying_axes[1]]) * static_cast<long>(AS)};
                const int face_local = node_index(AFaceGrid.local_nodes, ALocalNodes[c]);
                const auto &pq = FACE_CORNER_PQ[static_cast<std::size_t>(face_local)];
                block_to_face_pq[k] = {pq.first * static_cast<long>(AS), pq.second * static_cast<long>(AS)};
            }

            for (std::size_t i = 1; i < AS; ++i) {
                for (std::size_t j = 1; j < AS; ++j) {
                    const auto [u_idx, v_idx] = bilinear_index_remap(block_to_uvw, i, j, AS);
                    const auto [p, q] = bilinear_index_remap(block_to_face_pq, i, j, AS);
                    std::array<std::size_t, 3> idx{};
                    idx[fixed_axis] = fixed_value * AS;
                    idx[varying_axes[0]] = u_idx;
                    idx[varying_axes[1]] = v_idx;
                    AGrid[idx[0]][idx[1]][idx[2]] = AFaceGrid.grid[p][q];
                }
            }
        }

        /**
         * @brief Builds one hex block's `(S+1)^3` node-id grid (corners/edges/boundary-face-interiors
         * reused, only the `(S-1)^3` strictly interior points newly sampled from the block's own
         * volume) and emits its `S^3` sub-cells into the mesh.
         * @param ABlock The block to sample.
         * @param ANodeIds Already-populated node-attribute-to-mesh-id map.
         * @param AEdgeChains Already-populated edge-attribute-to-chain map.
         * @param AFaceGrids Already-populated face-attribute-to-grid map (its 6 boundary faces).
         * @param AMesh The mesh being built, appended to for interior points and sub-cells.
         * @param ADims The `NODE_CLASSIFICATION_DIM_VARIABLE` variable, extended for interior points.
         * @param ATags The `NODE_CLASSIFICATION_TAG_VARIABLE` variable, extended for interior points.
         * @param AS Subdivisions per axis.
         */
        void emit_hex_cells(Block ABlock,
                            const std::map<Node, NodeId> &ANodeIds,
                            const std::map<Edge, EdgeChain> &AEdgeChains,
                            const std::map<Face, FaceGrid> &AFaceGrids,
                            UnstructuredMesh<CubicTraits> &AMesh,
                            Variable<Int> &ADims,
                            Variable<Int> &ATags,
                            std::size_t AS) {
            const Dart bd = ABlock->dart();
            std::array<Node, 8> local_nodes{};
            local_nodes[0] = m_cmap.template attribute<0>(bd);
            local_nodes[1] = m_cmap.template attribute<0>(m_cmap.template beta<1>(bd));
            local_nodes[2] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1>(bd));
            local_nodes[3] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1, 1>(bd));
            local_nodes[4] = m_cmap.template attribute<0>(m_cmap.template beta<2, 1, 1>(bd));
            local_nodes[5] = m_cmap.template attribute<0>(m_cmap.template beta<1, 2, 1, 1>(bd));
            local_nodes[6] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1, 2, 1, 1>(bd));
            local_nodes[7] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1, 1, 2, 1, 1>(bd));

            Grid3D grid(AS + 1, std::vector<std::vector<NodeId>>(AS + 1, std::vector<NodeId>(AS + 1)));

            for (std::size_t c = 0; c < 8; ++c) {
                const auto &uvw = HEX_CORNER_UVW[c];
                grid[static_cast<std::size_t>(uvw[0]) * AS][static_cast<std::size_t>(uvw[1]) * AS]
                    [static_cast<std::size_t>(uvw[2]) * AS] = ANodeIds.at(local_nodes[c]);
            }

            for (auto it = m_cmap.template one_dart_per_incident_cell<1, 3>(bd).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<1, 3>(bd).end();
                 it != itend;
                 ++it) {
                const int a = node_index(local_nodes, m_cmap.template attribute<0>(it));
                const int b = node_index(local_nodes, m_cmap.template attribute<0>(m_cmap.template beta<1>(it)));
                const std::size_t k = find_edge(HEX_EDGES, a, b);
                const Edge e = m_cmap.template attribute<1>(it);
                fill_hex_edge_line(grid, HEX_EDGES[k], local_nodes, AEdgeChains.at(e), AS);
            }

            for (auto it = m_cmap.template one_dart_per_incident_cell<2, 3>(bd).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<2, 3>(bd).end();
                 it != itend;
                 ++it) {
                std::array<int, 4> corner_idx{};
                Dart wd = it;
                for (std::size_t c = 0; c < 4; ++c) {
                    corner_idx[c] = node_index(local_nodes, m_cmap.template attribute<0>(wd));
                    wd = m_cmap.template beta<1>(wd);
                }
                const Face face = m_cmap.template attribute<2>(it);
                fill_hex_face_interior(grid, corner_idx, local_nodes, AFaceGrids.at(face), AS);
            }

            for (std::size_t i = 1; i < AS; ++i) {
                for (std::size_t j = 1; j < AS; ++j) {
                    for (std::size_t k = 1; k < AS; ++k) {
                        const double u = static_cast<double>(i) / static_cast<double>(AS);
                        const double v = static_cast<double>(j) / static_cast<double>(AS);
                        const double w = static_cast<double>(k) / static_cast<double>(AS);
                        grid[i][j][k] = AMesh.add_node(ABlock->info().volume.value(u, v, w));
                        record_node_classification(ADims, ATags, grid[i][j][k], ABlock->info().geom_targets);
                    }
                }
            }

            for (std::size_t i = 0; i < AS; ++i) {
                for (std::size_t j = 0; j < AS; ++j) {
                    for (std::size_t k = 0; k < AS; ++k) {
                        AMesh.add_cell(grid[i][j][k],
                                       grid[i + 1][j][k],
                                       grid[i + 1][j + 1][k],
                                       grid[i][j + 1][k],
                                       grid[i][j][k + 1],
                                       grid[i + 1][j][k + 1],
                                       grid[i + 1][j + 1][k + 1],
                                       grid[i][j + 1][k + 1]);
                    }
                }
            }
        }

        /**
         * @brief Creates a node attribute at a given position.
         * @param APoint Spatial location of the node.
         * @return The created node.
         */
        Node create_node(const Point3d &APoint) {
            Node n = m_cmap.template create_attribute<0>();
            n->info().point = APoint;
            n->info().id = m_next_node_id++;
            return n;
        }

        /**
         * @brief Creates a default (unclassified) edge attribute; its curve is set separately.
         * @return The created edge.
         */
        Edge create_edge() { return m_cmap.template create_attribute<1>(); }

        /**
         * @brief Creates a default (unclassified) face attribute; its surface is set separately.
         * @return The created face.
         */
        Face create_face() { return m_cmap.template create_attribute<2>(); }

        /**
         * @brief Builds a straight `TEdgeCurve` between 2 points: `NumControlPoints` control
         * points, linearly interpolated.
         * @param AA Start point (parameter 0).
         * @param AB End point (parameter 1).
         * @return The straight curve.
         */
        TEdgeCurve straight_curve(const Point3d &AA, const Point3d &AB) const {
            TEdgeCurve curve(m_degree);
            const std::size_t n = m_degree + 1;
            const Vector3d ab(AA, AB);
            for (std::size_t i = 0; i < n; ++i) {
                const double t = (n > 1) ? static_cast<double>(i) / static_cast<double>(n - 1) : 0.0;
                curve[i] = AA + ab * t;
            }
            return curve;
        }

        /**
         * @brief Finds the index of the entry in @p AEdges whose (from,to) pair matches
         * `{a,b}`, in either order.
         * @tparam TN Number of entries in @p AEdges.
         * @param AEdges The edge table to search.
         * @param a First node index.
         * @param b Second node index.
         * @return The matching entry's index.
         * @pre A matching entry must exist.
         */
        template<std::size_t TN>
        static std::size_t find_edge(const std::array<std::pair<int, int>, TN> &AEdges, int a, int b) {
            for (std::size_t k = 0; k < TN; ++k) {
                if ((AEdges[k].first == a && AEdges[k].second == b) ||
                    (AEdges[k].first == b && AEdges[k].second == a)) {
                    return k;
                }
            }
            assert(false && "Blocking::find_edge: no matching edge found");
            return 0;
        }

        /**
         * @brief Finds the index of the entry in `HEX_FACES` whose corner set matches @p ACorners.
         * @param ACorners The 4 node indices bounding a face (any order).
         * @return The matching entry's index in `HEX_FACES`.
         * @pre A matching entry must exist.
         */
        static std::size_t find_face(const std::array<int, 4> &ACorners) {
            for (std::size_t f = 0; f < 6; ++f) {
                bool all_found = true;
                for (int c : ACorners) {
                    bool found = false;
                    for (int fc : HEX_FACES[f].corners) {
                        if (fc == c) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        all_found = false;
                        break;
                    }
                }
                if (all_found) return f;
            }
            assert(false && "Blocking::find_face: no matching face found");
            return 0;
        }

        /**
         * @brief Finds the index of @p AN in @p ANodes.
         * @param ANodes The node table to search.
         * @param AN The node to find.
         * @return The matching index.
         * @pre @p AN must be one of @p ANodes.
         */
        template<std::size_t TN>
        static int node_index(const std::array<Node, TN> &ANodes, Node AN) {
            for (std::size_t i = 0; i < TN; ++i) {
                if (ANodes[i] == AN) return static_cast<int>(i);
            }
            assert(false && "Blocking::node_index: node not found");
            return -1;
        }

        /**
         * @brief One entry per bounding face of a hex block: which 4 nodes it spans (for matching
         * a CGAL face-dart to this table, order-independent), and which 4 of `HEX_EDGES` form its
         * `EdgeU0`/`EdgeU1`/`EdgeV0`/`EdgeV1` boundary for `coons_surface_from_edges` (see
         * math/CoonsPatch.h's `[v][w]`/`[u][w]`/`[u][v]` indexing convention per face). Derived by
         * hand from `CubicTraits::Cell`'s own HEX8 corner layout and cross-checked end-to-end by
         * `blocking_creation_tests.cpp` (a built hex's TFI-evaluated center must be its geometric
         * centroid).
         */
        struct HexFaceSpec {
            /** @brief The 4 node indices bounding this face (unordered). */
            std::array<int, 4> corners;
            /** @brief Index into `HEX_EDGES` for this face's Coons `EdgeU0`. */
            std::size_t edge_u0;
            /** @brief Index into `HEX_EDGES` for this face's Coons `EdgeU1`. */
            std::size_t edge_u1;
            /** @brief Index into `HEX_EDGES` for this face's Coons `EdgeV0`. */
            std::size_t edge_v0;
            /** @brief Index into `HEX_EDGES` for this face's Coons `EdgeV1`. */
            std::size_t edge_v1;
        };

        /** @brief The 4 boundary edges of a quad face (standalone 2D block, or one of a hex's 6
         * bounding faces), as (from,to) local-node-index pairs, in the fixed direction
         * `coons_surface_from_edges()` needs: index 0 = `EdgeU0`, 2 = `EdgeU1`, 3 = `EdgeV0`, 1 =
         * `EdgeV1` (u along local corner 0->1, v along local corner 0->3). */
        static constexpr std::array<std::pair<int, int>, 4> QUAD_EDGES = {std::pair{0, 1},
                                                                          std::pair{1, 2},
                                                                          std::pair{3, 2},
                                                                          std::pair{0, 3}};

        /** @brief The 12 edges of a hex block, as (from,to) node-index pairs, in the single fixed
         * direction `HEX_FACES` always uses (verified consistent for both incident faces of every
         * edge — no reversal ever needed). Order: {0,3},{4,7},{0,4},{3,7},{1,2},{5,6},{1,5},{2,6},
         * {0,1},{4,5},{3,2},{7,6}. */
        static constexpr std::array<std::pair<int, int>, 12> HEX_EDGES = {std::pair{0, 3},
                                                                          std::pair{4, 7},
                                                                          std::pair{0, 4},
                                                                          std::pair{3, 7},
                                                                          std::pair{1, 2},
                                                                          std::pair{5, 6},
                                                                          std::pair{1, 5},
                                                                          std::pair{2, 6},
                                                                          std::pair{0, 1},
                                                                          std::pair{4, 5},
                                                                          std::pair{3, 2},
                                                                          std::pair{7, 6}};

        /** @brief The 6 bounding faces of a hex block: `{Fu0, Fu1, Fv0, Fv1, Fw0, Fw1}`, matching
         * `tfi_volume_from_faces()`'s expected argument order. */
        static constexpr std::array<HexFaceSpec, 6> HEX_FACES = {
            HexFaceSpec{{0, 3, 4, 7}, 0, 1, 2, 3},   // Fu0 (u=0): grid[v][w]
            HexFaceSpec{{1, 2, 5, 6}, 4, 5, 6, 7},   // Fu1 (u=1): grid[v][w]
            HexFaceSpec{{0, 1, 4, 5}, 8, 9, 2, 6},   // Fv0 (v=0): grid[u][w]
            HexFaceSpec{{3, 2, 7, 6}, 10, 11, 3, 7}, // Fv1 (v=1): grid[u][w]
            HexFaceSpec{{0, 1, 3, 2}, 8, 10, 0, 4},  // Fw0 (w=0): grid[u][v]
            HexFaceSpec{{4, 5, 7, 6}, 9, 11, 1, 5}   // Fw1 (w=1): grid[u][v]
        };

        // ---------------------------------------------------------------------
        // Sheet cut: collecting the sheet
        // ---------------------------------------------------------------------

        /**
         * @brief The key a node is ordered by when a canonical frame has to choose a starting corner.
         *
         * `NodeInfo::id` first, which survives CGAL rebuilding the attribute behind a node; then the
         * node's own coordinates, purely to break the tie that 2 attributes carrying the *copied* id
         * of the node they were split from would otherwise leave.
         *
         * That tie is not hypothetical — CGAL's split copies the info, id and all — and it must not
         * be broken on the attribute handle. A handle is an address, addresses move between runs,
         * and a frame resting on one answers differently every time the program starts. Coordinates
         * are at least the same on every run, and a tie surviving *them* is 2 nodes at one point,
         * where the frames they order are congruent anyway.
         *
         * @param ANode The node to order.
         * @return Its ordering key.
         */
        /** @brief Every face's canonical frame, as the positions of its 4 corners. */
        using FaceFrames = std::map<Face, std::array<Point3d, 4>>;
        /** @brief Every block's canonical frame, as the positions of its 8 corners. */
        using BlockFrames = std::map<Block, std::array<Point3d, 8>>;

        /**
         * @brief Records the canonical frame of every face and block, as corner *positions*.
         *
         * Taken before a topological edit so that `reframe_changed_cells()` can put back whatever the
         * edit moves. Positions rather than ids because the edit is exactly what changes ids: CGAL
         * rebuilds an attribute whose orbit came apart, and the node behind it is then a different
         * attribute needing a new id. Positions are the one thing a removal leaves alone — the values
         * compared are literally the same stored doubles, not 2 computations of one point — so a
         * cell's corners can still be recognised on the other side of it.
         *
         * @param AFaces Filled with each face's frame.
         * @param ABlocks Filled with each block's frame.
         */
        void record_frames(FaceFrames &AFaces, BlockFrames &ABlocks) {
            for (auto it = m_cmap.template attributes<2>().begin(), itend = m_cmap.template attributes<2>().end();
                 it != itend;
                 ++it) {
                const auto corners = face_local_nodes(it);
                std::array<Point3d, 4> frame{};
                for (std::size_t c = 0; c < 4; ++c) {
                    frame[c] = corners[c]->info().point;
                }
                AFaces.emplace(it, frame);
            }
            for (auto it = m_cmap.template attributes<3>().begin(), itend = m_cmap.template attributes<3>().end();
                 it != itend;
                 ++it) {
                const auto corners = block_local_nodes(it);
                std::array<Point3d, 8> frame{};
                for (std::size_t c = 0; c < 8; ++c) {
                    frame[c] = corners[c]->info().point;
                }
                ABlocks.emplace(it, frame);
            }
        }

        /**
         * @brief Re-expresses the stored geometry of every surviving cell whose canonical frame the
         * edit just moved.
         *
         * A frame is derived, not stored, so an edit that renumbers a corner silently rotates or
         * reflects it — and the surface or volume sitting in it then reads back transformed. Left
         * alone, the next cut through such a face takes its isocurve from the wrong side and leaves a
         * straight edge with a bow in it. Nothing is recomputed here: the geometry is permuted into
         * its cell's new frame, which is exact.
         *
         * @param AFaces Each face's frame before the edit, from `record_frames()`.
         * @param ABlocks Each block's frame before the edit.
         */
        void reframe_changed_cells(const FaceFrames &AFaces, const BlockFrames &ABlocks) {
            for (auto it = m_cmap.template attributes<2>().begin(), itend = m_cmap.template attributes<2>().end();
                 it != itend;
                 ++it) {
                const auto before = AFaces.find(it);
                if (before == AFaces.end()) continue;
                const auto corners = face_local_nodes(it);
                std::array<std::array<int, 2>, 4> corner_ab{};
                bool moved = false;
                bool reshaped = false;
                for (std::size_t c = 0; c < 4; ++c) {
                    const std::size_t was = position_index(before->second, corners[c]->info().point);
                    if (was == 4) {
                        // A corner the face did not have before: the edit split this face rather than
                        // just renumbering it, and whoever split it already stored its geometry in the
                        // frame it has now.
                        reshaped = true;
                        break;
                    }
                    corner_ab[c] = QUAD_CORNER_IJ[was];
                    if (was != c) moved = true;
                }
                if (!reshaped && moved) it->info().surface = reframed_surface(it->info().surface, corner_ab);
            }

            for (auto it = m_cmap.template attributes<3>().begin(), itend = m_cmap.template attributes<3>().end();
                 it != itend;
                 ++it) {
                const auto before = ABlocks.find(it);
                if (before == ABlocks.end()) continue;
                const auto corners = block_local_nodes(it);
                std::array<std::array<int, 3>, 8> corner_uvw{};
                bool moved = false;
                bool reshaped = false;
                for (std::size_t c = 0; c < 8; ++c) {
                    const std::size_t was = position_index(before->second, corners[c]->info().point);
                    if (was == 8) {
                        reshaped = true;
                        break;
                    }
                    corner_uvw[c] = HEX_CORNER_UVW[was];
                    if (was != c) moved = true;
                }
                if (!reshaped && moved) it->info().volume = reframed_volume(it->info().volume, corner_uvw);
            }
        }

        /**
         * @brief Where @p APoint sits in @p AFrame.
         * @tparam TN Number of corners.
         * @param AFrame A recorded frame's corner positions.
         * @param APoint The corner to locate.
         * @return Its index in @p AFrame, or `TN` when it is not one of its corners — which marks a
         *         cell the edit reshaped rather than merely renumbered.
         */
        template<std::size_t TN>
        static std::size_t position_index(const std::array<Point3d, TN> &AFrame, const Point3d &APoint) {
            for (std::size_t i = 0; i < TN; ++i) {
                if (AFrame[i] == APoint) return i;
            }
            return TN;
        }

        /**
         * @brief The key an edge is ordered by, as the ids of the 2 nodes it joins, smaller first.
         *
         * Gives the class a run-to-run stable order on edges, which it does not otherwise have —
         * edges carry no id of their own, and ordering them by attribute handle is ordering them by
         * address.
         *
         * @param AEdge The edge to order.
         * @return Its ordering key.
         */
        std::pair<Int, Int> edge_order(Edge AEdge) {
            const Dart d = AEdge->dart();
            const Int a = m_cmap.template attribute<0>(d)->info().id;
            const Int b = m_cmap.template attribute<0>(m_cmap.template beta<1>(d))->info().id;
            return {std::min(a, b), std::max(a, b)};
        }

        /**
         * @brief Runs one topological edit, giving ids to whatever nodes it creates and putting back
         * the frames it moved.
         *
         * The 3 things every edit here has to be wrapped in, in the one order that works: note the
         * frames, edit, hand out ids for the attributes CGAL rebuilt, then re-express the cells whose
         * frame those new ids moved. Skipping any of them leaves geometry sitting in a frame that is
         * no longer the one it will be read in.
         *
         * @tparam TEdit A callable performing the edit.
         * @param AEdit The edit to run.
         */
        template<typename TEdit>
        void run_reframed(TEdit AEdit) {
            FaceFrames faces_before;
            BlockFrames blocks_before;
            record_frames(faces_before, blocks_before);
            AEdit();
            assign_missing_node_ids();
            reframe_changed_cells(faces_before, blocks_before);
        }

        /**
         * @brief Gives an id to every node that has none.
         *
         * A node arrives without one when CGAL splits an existing attribute (see `NodeSplitFunctor`),
         * which happens inside the very cell insertions and removals this class is built on. Run at
         * the end of each of those, so no 2 nodes ever share an id and the order the canonical frames
         * rest on has no tie left to settle by luck.
         */
        void assign_missing_node_ids() {
            for (auto it = m_cmap.template attributes<0>().begin(), itend = m_cmap.template attributes<0>().end();
                 it != itend;
                 ++it) {
                if (it->info().id < 0) it->info().id = m_next_node_id++;
            }
        }

        static std::tuple<Int, double, double, double> node_order(Node ANode) {
            const Point3d &p = ANode->info().point;
            return {ANode->info().id, p.x(), p.y(), p.z()};
        }

        /**
         * @brief Records @p AEdge as part of a sheet, starting at @p AFrom.
         * @param AStart The edge-to-start-node map being built, modified in place.
         * @param APending The traversal work list, modified in place.
         * @param AEdge The edge to record.
         * @param AFrom The endpoint its parameter is measured from.
         * @return false if @p AEdge was already recorded starting at its *other* end — the sheet
         *         disagrees with itself and no single parameter cuts it homogeneously.
         */
        static bool extend_sheet(std::map<Edge, Node> &AStart, std::vector<Edge> &APending, Edge AEdge, Node AFrom) {
            const auto it = AStart.find(AEdge);
            if (it == AStart.end()) {
                AStart.emplace(AEdge, AFrom);
                APending.push_back(AEdge);
                return true;
            }
            return it->second == AFrom;
        }

        /** @brief The blocks incident to one edge.
         * @param AEdge The edge to inspect.
         * @return Its incident blocks (empty for an edge of a standalone quad). */
        std::vector<Block> blocks_of_edge(Edge AEdge) {
            std::vector<Block> blocks;
            const Dart ed = AEdge->dart();
            for (auto it = m_cmap.template one_dart_per_incident_cell<3, 1>(ed).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<3, 1>(ed).end();
                 it != itend;
                 ++it) {
                const Block b = m_cmap.template attribute<3>(it);
                if (b != nullptr) blocks.push_back(b);
            }
            return blocks;
        }

        /** @brief The faces incident to one edge that belong to no block — a standalone 2D quad
         * block's own face, a hex's bounding faces being reached through the block instead.
         * @param AEdge The edge to inspect.
         * @return Its incident block-less faces. */
        std::vector<Face> standalone_faces_of_edge(Edge AEdge) {
            std::vector<Face> faces;
            const Dart ed = AEdge->dart();
            for (auto it = m_cmap.template one_dart_per_incident_cell<2, 1>(ed).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<2, 1>(ed).end();
                 it != itend;
                 ++it) {
                if (m_cmap.template attribute<3>(it) != nullptr) continue;
                const Face f = m_cmap.template attribute<2>(it);
                if (f != nullptr) faces.push_back(f);
            }
            return faces;
        }

        /** @brief One block's 8 corners in `HEX_CORNER_UVW` order, walked from its own dart.
         * @param ABlock The block to inspect.
         * @return Its 8 corners.
         * @pre The block must still be a plain hex (no edge of it split yet). */
        std::array<Node, 8> block_local_nodes(Block ABlock) {
            // Read from whichever of the block's darts gives the lexicographically smallest tuple of
            // corner handles, never from `ABlock->dart()`. Every dart of a hex yields a valid frame —
            // the 24 of them are the cube's rotations — and which one the attribute names is CGAL's
            // to change: it re-points an attribute the moment the dart it named is erased, which a
            // deletion next door does. A frame that moves silently rotates everything stored in it,
            // and the block's volume is stored in it. See `face_local_nodes()` for the same trap one
            // dimension down, and what it does to a straight edge.
            std::array<Node, 8> best{};
            std::array<std::tuple<Int, double, double, double>, 8> best_key{};
            bool have = false;
            for (auto it = m_cmap.template darts_of_cell<3>(ABlock->dart()).begin(),
                      itend = m_cmap.template darts_of_cell<3>(ABlock->dart()).end();
                 it != itend;
                 ++it) {
                const Dart bd = it;
                std::array<Node, 8> n{};
                n[0] = m_cmap.template attribute<0>(bd);
                n[1] = m_cmap.template attribute<0>(m_cmap.template beta<1>(bd));
                n[2] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1>(bd));
                n[3] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1, 1>(bd));
                n[4] = m_cmap.template attribute<0>(m_cmap.template beta<2, 1, 1>(bd));
                n[5] = m_cmap.template attribute<0>(m_cmap.template beta<1, 2, 1, 1>(bd));
                n[6] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1, 2, 1, 1>(bd));
                n[7] = m_cmap.template attribute<0>(m_cmap.template beta<1, 1, 1, 2, 1, 1>(bd));
                std::array<std::tuple<Int, double, double, double>, 8> key{};
                for (std::size_t c = 0; c < 8; ++c) {
                    key[c] = node_order(n[c]);
                }
                if (!have || key < best_key) {
                    best = n;
                    best_key = key;
                    have = true;
                }
            }
            return best;
        }

        /** @brief One face's 4 corners in its own dart-cycle order.
         * @param AFace The face to inspect.
         * @return Its 4 corners.
         * @pre The face must still be a plain quad (no edge of it split yet). */
        std::array<Node, 4> face_local_nodes(Face AFace) {
            std::array<Node, 4> cycle{};
            Dart walk = AFace->dart();
            for (std::size_t c = 0; c < 4; ++c) {
                cycle[c] = m_cmap.template attribute<0>(walk);
                walk = m_cmap.template beta<1>(walk);
            }

            // Rotated to start at the lowest-handle corner, and run towards whichever of its 2
            // neighbours is lower — which pins the frame against both the rotation and the
            // reflection that `dart()` alone leaves free.
            //
            // Anchoring it to `dart()` instead is what bent straight edges after a deletion. A face
            // between 2 blocks carries darts on both sides, and the 2 sides walk its cycle in
            // *opposite* directions; deleting one block erases that side's darts, so CGAL re-points
            // the attribute to the surviving side and the frame silently mirrors. Everything stored
            // in the old frame — the surface above all — then reads back reflected, and the next cut
            // through that face takes its isocurve from the wrong side of the cut. Its 2 ends still
            // land on the right nodes, because they are pinned there, and its middle bulges out to
            // the mirrored position: a straight edge with a bow in it.
            std::size_t start = 0;
            for (std::size_t c = 1; c < 4; ++c) {
                if (node_order(cycle[c]) < node_order(cycle[start])) start = c;
            }
            const bool forward = node_order(cycle[(start + 1) % 4]) < node_order(cycle[(start + 3) % 4]);

            std::array<Node, 4> canonical{};
            for (std::size_t i = 0; i < 4; ++i) {
                canonical[i] = forward ? cycle[(start + i) % 4] : cycle[(start + 4 - i) % 4];
            }
            return canonical;
        }

        /** @brief One block's 12 edges, indexed to match `HEX_EDGES`.
         * @param ABlock The block to inspect.
         * @param ACorners Its 8 corners, from `block_local_nodes()`.
         * @return Its 12 edges in `HEX_EDGES` order. */
        std::array<Edge, 12> block_edges(Block ABlock, const std::array<Node, 8> &ACorners) {
            std::array<Edge, 12> edges{};
            const Dart bd = ABlock->dart();
            for (auto it = m_cmap.template one_dart_per_incident_cell<1, 3>(bd).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<1, 3>(bd).end();
                 it != itend;
                 ++it) {
                const int a = node_index(ACorners, m_cmap.template attribute<0>(it));
                const int b = node_index(ACorners, m_cmap.template attribute<0>(m_cmap.template beta<1>(it)));
                edges[find_edge(HEX_EDGES, a, b)] = m_cmap.template attribute<1>(it);
            }
            return edges;
        }

        /** @brief One face's 4 edges, indexed to match `QUAD_EDGES`.
         * @param AFace The face to inspect.
         * @param ACorners Its 4 corners, from `face_local_nodes()`.
         * @return Its 4 edges in `QUAD_EDGES` order. */
        std::array<Edge, 4> face_edges(Face AFace, const std::array<Node, 4> &ACorners) {
            std::array<Edge, 4> edges{};
            const Dart fd = AFace->dart();
            for (auto it = m_cmap.template one_dart_per_incident_cell<1, 2>(fd).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<1, 2>(fd).end();
                 it != itend;
                 ++it) {
                const int a = node_index(ACorners, m_cmap.template attribute<0>(it));
                const int b = node_index(ACorners, m_cmap.template attribute<0>(m_cmap.template beta<1>(it)));
                edges[find_edge(QUAD_EDGES, a, b)] = m_cmap.template attribute<1>(it);
            }
            return edges;
        }

        /** @brief Finds which slot of @p ATable holds @p AEdge.
         * @tparam TN Number of slots.
         * @param ATable The edge table to search.
         * @param AEdge The edge to find.
         * @return Its index.
         * @pre @p AEdge must be in @p ATable. */
        template<std::size_t TN>
        static std::size_t table_index(const std::array<Edge, TN> &ATable, Edge AEdge) {
            for (std::size_t k = 0; k < TN; ++k) {
                if (ATable[k] == AEdge) return k;
            }
            assert(false && "Blocking::table_index: edge not found");
            return 0;
        }

        /** @brief Which local axis one of `HEX_EDGES`' 12 entries runs along.
         * @param AK Index into `HEX_EDGES`.
         * @return 0/1/2 for `u`/`v`/`w`. */
        static int hex_edge_axis(std::size_t AK) {
            const auto a = static_cast<std::size_t>(HEX_EDGES[AK].first);
            const auto b = static_cast<std::size_t>(HEX_EDGES[AK].second);
            for (int ax = 0; ax < 3; ++ax) {
                if (HEX_CORNER_UVW[a][static_cast<std::size_t>(ax)] !=
                    HEX_CORNER_UVW[b][static_cast<std::size_t>(ax)]) {
                    return ax;
                }
            }
            assert(false && "Blocking::hex_edge_axis: degenerate edge");
            return 0;
        }

        /** @brief Which local axis one of `QUAD_EDGES`' 4 entries runs along.
         * @param AK Index into `QUAD_EDGES`.
         * @return 0 for `u`, 1 for `v`. */
        static int quad_edge_axis(std::size_t AK) {
            const auto a = static_cast<std::size_t>(QUAD_EDGES[AK].first);
            const auto b = static_cast<std::size_t>(QUAD_EDGES[AK].second);
            return (QUAD_CORNER_IJ[a][0] != QUAD_CORNER_IJ[b][0]) ? 0 : 1;
        }

        /**
         * @brief Adds to @p ASheet every face crossed by the edges of @p AStart, rejecting any face
         * the sheet meets in a way that leaves no single cut.
         * @param AStart The sheet's edge-to-start-node map.
         * @param ASheet The sheet being built, modified in place.
         * @return false if some face carries a number of sheet edges other than 0 or 2, or carries 2
         *         that are not parallel to each other.
         */
        bool collect_sheet_faces(const std::map<Edge, Node> &AStart, Sheet &ASheet) {
            for (auto it = m_cmap.template attributes<2>().begin(), itend = m_cmap.template attributes<2>().end();
                 it != itend;
                 ++it) {
                const Face f = it;
                const std::array<Node, 4> corners = face_local_nodes(f);
                const std::array<Edge, 4> edges = face_edges(f, corners);

                std::vector<std::size_t> hits;
                for (std::size_t k = 0; k < 4; ++k) {
                    if (AStart.count(edges[k]) != 0) hits.push_back(k);
                }
                if (hits.empty()) continue;
                if (hits.size() != 2) return false;
                const int axis = quad_edge_axis(hits[0]);
                if (quad_edge_axis(hits[1]) != axis) return false;

                const Node from = AStart.at(edges[hits[0]]);
                ASheet.faces.push_back(SheetFace{f, corners, axis, node_index(corners, from)});
            }
            return true;
        }

        /**
         * @brief Adds to @p ASheet every block crossed by the edges of @p AStart, rejecting any block
         * the sheet meets in a way that leaves no single cut.
         * @param AStart The sheet's edge-to-start-node map.
         * @param ASheet The sheet being built, modified in place.
         * @return false if some block carries a number of sheet edges other than 0 or 4, or carries 4
         *         that are not all parallel to each other.
         */
        bool collect_sheet_blocks(const std::map<Edge, Node> &AStart, Sheet &ASheet) {
            for (auto it = m_cmap.template attributes<3>().begin(), itend = m_cmap.template attributes<3>().end();
                 it != itend;
                 ++it) {
                const Block b = it;
                const std::array<Node, 8> corners = block_local_nodes(b);
                const std::array<Edge, 12> edges = block_edges(b, corners);

                std::vector<std::size_t> hits;
                for (std::size_t k = 0; k < 12; ++k) {
                    if (AStart.count(edges[k]) != 0) hits.push_back(k);
                }
                if (hits.empty()) continue;
                if (hits.size() != 4) return false;
                const int axis = hex_edge_axis(hits[0]);
                for (const std::size_t k : hits) {
                    if (hex_edge_axis(k) != axis) return false;
                }

                const Node from = AStart.at(edges[hits[0]]);
                ASheet.blocks.push_back(SheetBlock{b, corners, axis, node_index(corners, from)});
            }
            return true;
        }

        // ---------------------------------------------------------------------
        // Sheet cut: applying the cut
        // ---------------------------------------------------------------------

        /** @brief A new node inserted mid-edge by a cut, and the 2 parent corners it sits between. */
        struct MidNode {
            /** @brief The parent endpoint the cut parameter was measured from. */
            Node from{};
            /** @brief The parent endpoint at the other end. */
            Node to{};
        };

        /** @brief Like `node_index()`, but reports absence instead of asserting.
         * @tparam TN Number of slots.
         * @param ANodes The node table to search.
         * @param AN The node to find.
         * @return Its index, or -1 when @p AN is not in @p ANodes. */
        template<std::size_t TN>
        static int node_index_or_none(const std::array<Node, TN> &ANodes, Node AN) {
            for (std::size_t i = 0; i < TN; ++i) {
                if (ANodes[i] == AN) return static_cast<int>(i);
            }
            return -1;
        }

        /** @brief A dart of @p AEdge whose source node is @p AFrom.
         * @param AEdge The edge to look in.
         * @param AFrom The node the dart must start at.
         * @return Such a dart.
         * @pre @p AFrom must be one of @p AEdge's 2 endpoints. */
        Dart dart_of_edge_from(Edge AEdge, Node AFrom) {
            for (auto it = m_cmap.template darts_of_cell<1>(AEdge->dart()).begin(),
                      itend = m_cmap.template darts_of_cell<1>(AEdge->dart()).end();
                 it != itend;
                 ++it) {
                if (m_cmap.template attribute<0>(it) == AFrom) return it;
            }
            assert(false && "Blocking::dart_of_edge_from: node is not an endpoint of the edge");
            return AEdge->dart();
        }

        /**
         * @brief Moves a curve's 2 end control points exactly onto the 2 nodes it is meant to join.
         *
         * The curve an edge stores has to pass through that edge's own endpoint nodes exactly, and a
         * curve read off a surface's control grid only gets there to rounding: the same corner
         * reached through a surface and through a curve are equal as real numbers and not as doubles.
         * Left alone, that gap is invisible geometrically and lethal to the orientation tests this
         * class runs by comparing positions (`oriented_curve()`, `curve_start_node()`), which then
         * pick the wrong end of the edge.
         *
         * Which end is which is the caller's to know — it comes out of the frame the curve was read
         * from, not out of comparing coordinates.
         *
         * @param ACurve The curve to pin, already running from @p AStart towards @p AEnd.
         * @param AStart The position its first control point must take.
         * @param AEnd The position its last control point must take.
         * @return @p ACurve with its 2 end control points replaced by those positions.
         */
        static TEdgeCurve pinned_curve(TEdgeCurve ACurve, const Point3d &AStart, const Point3d &AEnd) {
            const std::size_t last = ACurve.nb_control_points() - 1;
            ACurve[0] = AStart;
            ACurve[last] = AEnd;
            return ACurve;
        }

        /**
         * @brief The isoparametric boundary curve of a surface at the far end of one axis — for the
         * low half of a split, that is exactly the curve the cut ran along.
         * @param ASurface The surface to read the boundary off.
         * @param AAxis 0 for the `u=1` boundary (parameterized by `v`), 1 for `v=1` (by `u`).
         * @return The boundary curve, exactly (it is a row/column of the control grid).
         */
        static TEdgeCurve surface_boundary_curve(const FaceSurfaceT &ASurface, int AAxis) {
            const std::size_t n = ASurface.degree();
            std::vector<Point3d> points(n + 1);
            for (std::size_t a = 0; a <= n; ++a) {
                points[a] = (AAxis == 0) ? ASurface.control_points()[n][a] : ASurface.control_points()[a][n];
            }
            return TEdgeCurve(points);
        }

        /**
         * @brief The isoparametric boundary surface of a volume at the far end of one axis — for the
         * low half of a split, exactly the surface the cut ran along.
         * @param AVolume The volume to read the boundary off.
         * @param AAxis 0/1/2 for the `u=1`/`v=1`/`w=1` boundary.
         * @return The boundary surface, exactly (it is a slice of the control grid), parameterized by
         *         the 2 remaining axes in increasing order.
         */
        static FaceSurfaceT volume_boundary_surface(const BlockVolumeT &AVolume, int AAxis) {
            const std::size_t n = AVolume.degree();
            auto grid = FaceSurfaceT::make_grid(n);
            for (std::size_t a = 0; a <= n; ++a) {
                for (std::size_t b = 0; b <= n; ++b) {
                    if (AAxis == 0) {
                        grid[a][b] = AVolume.control_points()[n][a][b];
                    } else if (AAxis == 1) {
                        grid[a][b] = AVolume.control_points()[a][n][b];
                    } else {
                        grid[a][b] = AVolume.control_points()[a][b][n];
                    }
                }
            }
            return FaceSurfaceT(grid);
        }

        /**
         * @brief Re-expresses a surface given in some source frame into a target frame that differs
         * from it by one of the 8 square symmetries.
         *
         * A cut produces sub-surfaces in the frame of whatever they were cut out of, but the face
         * attribute they land on has its own frame, fixed by its dart cycle and not by us — and
         * everything downstream (`to_mesh()`, `rebuild_face_surface()`) reads a face's surface in that
         * own frame. Since both frames have the same 4 corners, the permutation is read straight off
         * where those corners sit, with no geometry involved.
         *
         * @param ASource The surface, in the source frame.
         * @param ACornerAB Where the target frame's local corners 0..3 sit in the source frame, each
         *        coordinate 0 or 1, in the target's own dart-cycle order.
         * @return The same surface, re-indexed into the target frame.
         */
        static FaceSurfaceT reframed_surface(const FaceSurfaceT &ASource,
                                             const std::array<std::array<int, 2>, 4> &ACornerAB) {
            const std::size_t n = ASource.degree();
            // Corner 0 -> corner 1 is the target's own u; whichever source axis moves along it is the
            // one the target's u maps to, and the remaining source axis takes the target's v.
            const bool a_follows_u = (ACornerAB[1][0] != ACornerAB[0][0]);
            auto grid = FaceSurfaceT::make_grid(n);
            for (std::size_t p = 0; p <= n; ++p) {
                for (std::size_t q = 0; q <= n; ++q) {
                    const std::size_t along_a = a_follows_u ? p : q;
                    const std::size_t along_b = a_follows_u ? q : p;
                    const std::size_t r = (ACornerAB[0][0] == 0) ? along_a : (n - along_a);
                    const std::size_t s = (ACornerAB[0][1] == 0) ? along_b : (n - along_b);
                    grid[p][q] = ASource.control_points()[r][s];
                }
            }
            return FaceSurfaceT(grid);
        }

        /**
         * @brief Splits every edge of a sheet, inserting one node per edge at the cut.
         * @param ASheet The sheet being cut.
         * @param AParam The cut parameter, as `cut_sheet()` defines it.
         * @param AMids Filled with every inserted node, mapped to the parent corners it lies between.
         */
        void split_sheet_edges(const Sheet &ASheet, double AParam, std::map<Node, MidNode> &AMids) {
            for (const SheetEdge &se : ASheet.edges) {
                const bool curve_starts_at_from = (curve_start_node(se.edge) == se.from);
                const double s = cut_parameter(se, AParam);
                const TEdgeCurve parent = se.edge->info().curve;
                const auto [low, high] = parent.split(s);

                Node mid = create_node(parent.value(s));
                // A node born on an edge lies on whatever that edge lies on.
                mid->info().geom_targets = se.edge->info().geom_targets;

                const Dart d = dart_of_edge_from(se.edge, se.from);
                const Node other = m_cmap.template attribute<0>(m_cmap.template beta<1>(d));
                const Dart nd = m_cmap.insert_cell_0_in_cell_1(d, mid);

                // insert_cell_0_in_cell_1 leaves `d` on the [from, mid] half and returns beta<1>(d),
                // which carries the [mid, other] one; both inherit the parent's classification
                // through CellData.h's SplitFunctor, and only the geometry needs setting. Which half
                // is which follows from where the parent's curve started, a node identity read before
                // the split — never from where the halves happen to land.
                store_curve(m_cmap.template attribute<1>(d),
                            curve_starts_at_from ? low : high,
                            curve_starts_at_from ? se.from : mid);
                store_curve(m_cmap.template attribute<1>(nd),
                            curve_starts_at_from ? high : low,
                            curve_starts_at_from ? mid : other);

                AMids.emplace(mid, MidNode{se.from, other});
            }
        }

        /**
         * @brief Splits every face of a sheet, inserting one edge per face along the cut.
         * @param ASheet The sheet being cut.
         * @param AParam The cut parameter, as `cut_sheet()` defines it.
         * @param AMids Every node `split_sheet_edges()` inserted.
         */
        void split_sheet_faces(const Sheet &ASheet, double AParam, const std::map<Node, MidNode> &AMids) {
            for (const SheetFace &sf : ASheet.faces) {
                const auto axis = static_cast<std::size_t>(sf.axis);
                const std::size_t other = 1 - axis;
                // The sheet's parameter is measured from `from_corner`; read in the face's own frame
                // it runs backwards whenever that corner sits at the far end of the cut axis.
                const bool forward = (QUAD_CORNER_IJ[static_cast<std::size_t>(sf.from_corner)][axis] == 0);
                const double t = forward ? AParam : 1.0 - AParam;

                const FaceSurfaceT parent = sf.face->info().surface;
                const auto [low, high] = (axis == 0) ? parent.split_u(t) : parent.split_v(t);

                // The 2 darts the new edge must span: the ones starting at this face's 2 new nodes.
                std::vector<Dart> mid_darts;
                Dart walk = sf.face->dart();
                const Dart first = walk;
                do {
                    if (AMids.count(m_cmap.template attribute<0>(walk)) != 0) mid_darts.push_back(walk);
                    walk = m_cmap.template beta<1>(walk);
                } while (walk != first);
                assert(mid_darts.size() == 2 && "Blocking::split_sheet_faces: face must carry exactly 2 cut nodes");

                const Dart nd = m_cmap.insert_cell_1_in_cell_2(mid_darts[0], mid_darts[1]);

                Edge new_edge = m_cmap.template attribute<1>(nd);
                if (new_edge == nullptr) {
                    new_edge = create_edge();
                    m_cmap.template set_attribute<1>(nd, new_edge);
                }
                // The cut runs along the low half's far isoparametric curve, exactly.
                // surface_boundary_curve() always runs along the axis the cut left alone, from its 0
                // end to its 1 end, so which of the 2 cut nodes it starts at is known outright — the
                // one born on the edge whose parent corners sit at 0 along that axis.
                const Node mid_a = m_cmap.template attribute<0>(mid_darts[0]);
                const Node mid_b = m_cmap.template attribute<0>(mid_darts[1]);
                const int side_a =
                    QUAD_CORNER_IJ[static_cast<std::size_t>(node_index(sf.corners, AMids.at(mid_a).from))][other];
                const Node curve_start = (side_a == 0) ? mid_a : mid_b;
                const Node curve_end = (side_a == 0) ? mid_b : mid_a;
                // Pinned onto those 2 nodes: read off a surface, the curve reaches them only to
                // rounding, and an edge's curve has to pass through its own endpoints exactly.
                store_curve(new_edge,
                            pinned_curve(surface_boundary_curve(low, sf.axis),
                                         curve_start->info().point,
                                         curve_end->info().point),
                            curve_start);
                // An edge born inside a face lies on whatever that face lies on.
                new_edge->info().geom_targets = sf.face->info().geom_targets;

                const Face half_a = m_cmap.template attribute<2>(mid_darts[0]);
                const Face half_b = m_cmap.template attribute<2>(mid_darts[1]);
                // Whichever half kept a parent corner from the near side of the cut is the low one.
                std::size_t low_corner = 0;
                for (std::size_t c = 0; c < 4; ++c) {
                    if (QUAD_CORNER_IJ[c][axis] == 0) low_corner = c;
                }
                const bool a_is_low = face_has_node(half_a, sf.corners[low_corner]);
                const Face low_face = a_is_low ? half_a : half_b;
                const Face high_face = a_is_low ? half_b : half_a;

                assign_half_face(low_face, low, sf, AMids, true, axis, other);
                assign_half_face(high_face, high, sf, AMids, false, axis, other);
            }
        }

        /**
         * @brief Stores one half of a split face's surface, re-expressed in that half's own frame.
         * @param AHalf The half face to write to.
         * @param ASurface Its surface, still in the parent face's frame.
         * @param ASheetFace The parent face, as the sheet captured it.
         * @param AMids Every node `split_sheet_edges()` inserted.
         * @param ALow Whether @p AHalf is the near-side half of the cut.
         * @param AAxis The parent's axis the cut ran across.
         * @param AOther The parent's other axis.
         */
        void assign_half_face(Face AHalf,
                              const FaceSurfaceT &ASurface,
                              const SheetFace &ASheetFace,
                              const std::map<Node, MidNode> &AMids,
                              bool ALow,
                              std::size_t AAxis,
                              std::size_t AOther) {
            const std::array<Node, 4> corners = face_local_nodes(AHalf);
            std::array<std::array<int, 2>, 4> corner_ab{};
            for (std::size_t c = 0; c < 4; ++c) {
                const auto mid = AMids.find(corners[c]);
                if (mid != AMids.end()) {
                    // A cut node sits at the far end of the near half and the near end of the far
                    // one; across the cut it keeps the parent coordinate of the edge it was born on.
                    const int idx = node_index(ASheetFace.corners, mid->second.from);
                    corner_ab[c][AAxis] = ALow ? 1 : 0;
                    corner_ab[c][AOther] = QUAD_CORNER_IJ[static_cast<std::size_t>(idx)][AOther];
                } else {
                    const int idx = node_index(ASheetFace.corners, corners[c]);
                    corner_ab[c] = QUAD_CORNER_IJ[static_cast<std::size_t>(idx)];
                }
            }
            AHalf->info().surface = reframed_surface(ASurface, corner_ab);
        }

        /**
         * @brief Re-expresses a volume given in some source frame into a target frame that differs
         * from it by one of the 48 cube symmetries.
         *
         * The 3D counterpart of `reframed_surface()`, and needed for the same reason: a cut produces
         * sub-volumes in the frame of the block they were cut out of, but each half is a *new* 3-cell
         * whose own frame comes from its own dart cycle — which CGAL builds, not us, and which is
         * generally some rotation or reflection of the parent's. Everything that later reads a
         * block's volume (`emit_hex_cells()`, `rebuild_block_volume()`) anchors `(u,v,w)` to that own
         * cycle, so the sub-volume has to be re-indexed into it.
         *
         * @param ASource The volume, in the source frame.
         * @param ACornerUVW Where the target frame's 8 local corners sit in the source frame, each
         *        coordinate 0 or 1, in `HEX_CORNER_UVW` order.
         * @return The same volume, re-indexed into the target frame.
         */
        static BlockVolumeT reframed_volume(const BlockVolumeT &ASource,
                                            const std::array<std::array<int, 3>, 8> &ACornerUVW) {
            const std::size_t n = ASource.degree();
            // Corner 0 to corners 1/3/4 walk the target's own u/v/w; whichever source axis moves
            // along each is the one it maps to, and corner 0's own position gives the direction.
            std::array<std::size_t, 3> axis_of{};
            for (std::size_t target = 0; target < 3; ++target) {
                const std::size_t neighbour = (target == 0) ? 1 : (target == 1 ? 3 : 4);
                for (std::size_t k = 0; k < 3; ++k) {
                    if (ACornerUVW[neighbour][k] != ACornerUVW[0][k]) axis_of[target] = k;
                }
            }

            auto grid = BlockVolumeT::make_grid(n);
            for (std::size_t p = 0; p <= n; ++p) {
                for (std::size_t q = 0; q <= n; ++q) {
                    for (std::size_t r = 0; r <= n; ++r) {
                        const std::array<std::size_t, 3> along{p, q, r};
                        std::array<std::size_t, 3> src{};
                        for (std::size_t target = 0; target < 3; ++target) {
                            const std::size_t k = axis_of[target];
                            src[k] = (ACornerUVW[0][k] == 0) ? along[target] : (n - along[target]);
                        }
                        grid[p][q][r] = ASource.control_points()[src[0]][src[1]][src[2]];
                    }
                }
            }
            return BlockVolumeT(grid);
        }

        /**
         * @brief Stores one half of a split block's volume, re-expressed in that half's own frame.
         * @param AHalf The half block to write to.
         * @param AVolume Its volume, still in the parent block's frame.
         * @param ASheetBlock The parent block, as the sheet captured it.
         * @param AMids Every node `split_sheet_edges()` inserted.
         * @param ALow Whether @p AHalf is the near-side half of the cut.
         * @param AAxis The parent's axis the cut ran across.
         */
        void assign_half_block(Block AHalf,
                               const BlockVolumeT &AVolume,
                               const SheetBlock &ASheetBlock,
                               const std::map<Node, MidNode> &AMids,
                               bool ALow,
                               std::size_t AAxis) {
            const std::array<Node, 8> corners = block_local_nodes(AHalf);
            std::array<std::array<int, 3>, 8> corner_uvw{};
            for (std::size_t c = 0; c < 8; ++c) {
                const auto mid = AMids.find(corners[c]);
                if (mid != AMids.end()) {
                    // A cut node keeps the parent coordinates of the edge it was born on across the
                    // 2 axes the cut left alone, and sits at the far end of the near half along the
                    // third — the near end of the far half.
                    const int idx = node_index(ASheetBlock.corners, mid->second.from);
                    corner_uvw[c] = HEX_CORNER_UVW[static_cast<std::size_t>(idx)];
                    corner_uvw[c][AAxis] = ALow ? 1 : 0;
                } else {
                    const int idx = node_index(ASheetBlock.corners, corners[c]);
                    corner_uvw[c] = HEX_CORNER_UVW[static_cast<std::size_t>(idx)];
                }
            }
            AHalf->info().volume = reframed_volume(AVolume, corner_uvw);
        }

        /**
         * @brief Splits every block of a sheet, inserting one face per block along the cut.
         * @param ASheet The sheet being cut.
         * @param AParam The cut parameter, as `cut_sheet()` defines it.
         * @param AMids Every node `split_sheet_edges()` inserted.
         */
        void split_sheet_blocks(const Sheet &ASheet, double AParam, const std::map<Node, MidNode> &AMids) {
            for (const SheetBlock &sb : ASheet.blocks) {
                const auto axis = static_cast<std::size_t>(sb.axis);
                const bool forward = (HEX_CORNER_UVW[static_cast<std::size_t>(sb.from_corner)][axis] == 0);
                const double t = forward ? AParam : 1.0 - AParam;

                const BlockVolumeT parent = sb.block->info().volume;
                const auto [low, high] = (axis == 0)   ? parent.split_u(t)
                                         : (axis == 1) ? parent.split_v(t)
                                                       : parent.split_w(t);

                // This block's own 4 cut nodes: those lying between 2 of its corners, along its cut
                // axis. Every other new node in the map belongs to some other block of the sheet.
                std::vector<Node> mids;
                for (const auto &[mid, ends] : AMids) {
                    const int ia = node_index_or_none(sb.corners, ends.from);
                    const int ib = node_index_or_none(sb.corners, ends.to);
                    if (ia < 0 || ib < 0) continue;
                    if (HEX_CORNER_UVW[static_cast<std::size_t>(ia)][axis] ==
                        HEX_CORNER_UVW[static_cast<std::size_t>(ib)][axis]) {
                        continue;
                    }
                    mids.push_back(mid);
                }
                // Sorted for the same reason the sheet's edges are: `AMids` is keyed on handles, and
                // `cut_loop()` starts from whichever of these comes first.
                std::sort(mids.begin(), mids.end(), [](Node a, Node b) { return a->info().id < b->info().id; });
                assert(mids.size() == 4 && "Blocking::split_sheet_blocks: block must carry exactly 4 cut nodes");

                const std::vector<Dart> path = cut_loop(sb.block, mids);
                const Dart nd = m_cmap.insert_cell_2_in_cell_3(path.begin(), path.end());

                Face new_face = m_cmap.template attribute<2>(nd);
                if (new_face == nullptr) {
                    new_face = create_face();
                    m_cmap.template set_attribute<2>(nd, new_face);
                }
                // A face born inside a block lies on whatever that block lies on.
                new_face->info().geom_targets = sb.block->info().geom_targets;
                assign_cut_face(new_face, volume_boundary_surface(low, sb.axis), sb, AMids, axis);

                const Block half_a = m_cmap.template attribute<3>(nd);
                const Block half_b = m_cmap.template attribute<3>(m_cmap.template beta<3>(nd));
                std::size_t low_corner = 0;
                for (std::size_t c = 0; c < 8; ++c) {
                    if (HEX_CORNER_UVW[c][axis] == 0) low_corner = c;
                }
                const bool a_is_low = block_has_node(half_a, sb.corners[low_corner]);
                assign_half_block(a_is_low ? half_a : half_b, low, sb, AMids, true, axis);
                assign_half_block(a_is_low ? half_b : half_a, high, sb, AMids, false, axis);
            }
        }

        /**
         * @brief Walks the closed 4-edge loop the cut traces around one block, as the path
         * `insert_cell_2_in_cell_3()` needs.
         * @param ABlock The block being split.
         * @param AMids Its own 4 cut nodes.
         * @return 4 darts, each running from one cut node to the next around the loop.
         */
        std::vector<Dart> cut_loop(Block ABlock, const std::vector<Node> &AMids) {
            std::vector<Dart> path;
            Node current = AMids[0];
            Node previous{};
            for (std::size_t step = 0; step < AMids.size(); ++step) {
                bool found = false;
                for (auto it = m_cmap.template darts_of_cell<3>(ABlock->dart()).begin(),
                          itend = m_cmap.template darts_of_cell<3>(ABlock->dart()).end();
                     it != itend && !found;
                     ++it) {
                    if (m_cmap.template attribute<0>(it) != current) continue;
                    const Node target = m_cmap.template attribute<0>(m_cmap.template beta<1>(it));
                    if (std::find(AMids.begin(), AMids.end(), target) == AMids.end()) continue;
                    if (step > 0 && target == previous) continue;
                    path.push_back(it);
                    previous = current;
                    current = target;
                    found = true;
                }
                assert(found && "Blocking::cut_loop: the cut nodes do not form a closed loop");
            }
            assert(current == AMids[0] && "Blocking::cut_loop: the loop did not close");
            return path;
        }

        /**
         * @brief Stores the surface of the face a block cut created, re-expressed in that face's own
         * frame.
         * @param AFace The new face.
         * @param ASurface Its surface, in the parent block's 2 non-cut axes, in increasing order.
         * @param ASheetBlock The parent block, as the sheet captured it.
         * @param AMids Every node `split_sheet_edges()` inserted.
         * @param AAxis The parent block's axis the cut ran across.
         */
        void assign_cut_face(Face AFace,
                             const FaceSurfaceT &ASurface,
                             const SheetBlock &ASheetBlock,
                             const std::map<Node, MidNode> &AMids,
                             std::size_t AAxis) {
            // volume_boundary_surface() parameterizes by the 2 axes the cut left alone, in order.
            std::array<std::size_t, 2> kept{};
            std::size_t k = 0;
            for (std::size_t ax = 0; ax < 3; ++ax) {
                if (ax != AAxis) kept[k++] = ax;
            }

            const std::array<Node, 4> corners = face_local_nodes(AFace);
            std::array<std::array<int, 2>, 4> corner_ab{};
            for (std::size_t c = 0; c < 4; ++c) {
                // Every corner of this face is a cut node, and keeps the parent-block coordinates of
                // the edge it was born on along both axes the cut left alone.
                const int idx = node_index(ASheetBlock.corners, AMids.at(corners[c]).from);
                for (std::size_t a = 0; a < 2; ++a) {
                    corner_ab[c][a] = HEX_CORNER_UVW[static_cast<std::size_t>(idx)][kept[a]];
                }
            }
            AFace->info().surface = reframed_surface(ASurface, corner_ab);
        }

        /**
         * @brief Rebuilds `m_hex_faces` from the blocks that currently exist.
         *
         * `m_hex_faces` means "every face belonging to some block": those are the candidates
         * `build_connectivity()` may sew at dimension 3, and the ones `to_mesh()` must *not* emit a
         * quad for, a quad being what a standalone 2D block is made of. Derived rather than
         * maintained by hand, because the operations that change it — a cut splitting a face in 2, a
         * deletion taking faces away with the block that owned them — would otherwise each have to
         * remember to, and a deletion cannot: once CGAL has garbage-collected a face attribute there
         * is no handle left to look up and drop.
         */
        void refresh_hex_faces() {
            std::set<Face> faces;
            for (auto it = m_cmap.template attributes<3>().begin(), itend = m_cmap.template attributes<3>().end();
                 it != itend;
                 ++it) {
                for (auto fit = m_cmap.template one_dart_per_incident_cell<2, 3>(it->dart()).begin(),
                          fend = m_cmap.template one_dart_per_incident_cell<2, 3>(it->dart()).end();
                     fit != fend;
                     ++fit) {
                    faces.insert(m_cmap.template attribute<2>(fit));
                }
            }
            m_hex_faces.assign(faces.begin(), faces.end());
        }

        /** @brief The degree every cell's geometry is built at — see `set_degree()`. */
        std::size_t m_degree = 1;
        /** @brief Source of `NodeInfo::id`. Only ever increases, so an id is never reused and the
         * order it defines never changes meaning. */
        Int m_next_node_id = 0;
        /** @brief Non-owning pointer to the geometric model this blocking is built against. */
        const TGeomModel *m_geom_model;
        /** @brief The underlying (always dimension-3) combinatorial map. */
        Map m_cmap;
        /** @brief Every face created as one of a hex block's 6 boundary faces (dim-3 sew candidates). */
        std::vector<Face> m_hex_faces;
    };

} // namespace gecko
