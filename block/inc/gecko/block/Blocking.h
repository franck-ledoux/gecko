#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <map>
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
     * `TEdgeCurve` (default `BezierCurve<1, Point3d>`, i.e. straight edges) is the "linear vs
     * curved blocking via template instantiation" axis from issue #22: any degree `N` gives curved
     * blocking, `N=1` collapses exactly to linear blocking through the same Coons/TFI construction
     * (see math/CoonsPatch.h). `TGeomModel` (`GeomModelConcept`-conforming, e.g. `FacetedGeometry`)
     * is the geometric model this blocking is built against, held for the whole lifetime of the
     * object — "a block structure always exists along a geometric model".
     *
     * @tparam TGeomModel Geometric model type, must satisfy `GeomModelConcept`.
     * @tparam TEdgeCurve Edge curve representation, must satisfy `EdgeCurveConcept`.
     */
    template<GeomModelConcept TGeomModel, EdgeCurveConcept TEdgeCurve = BezierCurve<1, Point3d>>
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
         * @brief Constructor.
         * @param AGeomModel The geometric model this blocking is built against; must outlive this
         *        Blocking (only a non-owning pointer is stored).
         */
        explicit Blocking(const TGeomModel &AGeomModel) : m_geom_model(&AGeomModel) {}

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
                e->info().curve = curves[k];
                m_cmap.template set_attribute<1>(it, e);
            }

            Face f = m_cmap.template create_attribute<2>();
            // EdgeU0, EdgeU1, EdgeV0, EdgeV1 per QUAD_EDGES' own indices 0,2,3,1.
            f->info().surface = coons_surface_from_edges(curves[0], curves[2], curves[3], curves[1]);
            m_cmap.template set_attribute<2>(d1, f);
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
                e->info().curve = curves[k];
                m_cmap.template set_attribute<1>(it, e);
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
                Face face = create_face();
                face->info().surface = face_surfaces[f];
                m_cmap.template set_attribute<2>(it, face);
                m_hex_faces.push_back(face);
            }

            Block b = m_cmap.template create_attribute<3>();
            b->info().volume = tfi_volume_from_faces(face_surfaces[0],
                                                     face_surfaces[1],
                                                     face_surfaces[2],
                                                     face_surfaces[3],
                                                     face_surfaces[4],
                                                     face_surfaces[5]);
            m_cmap.template set_attribute<3>(d1, b);
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
         * For each cell of topological dimension `d`, searches `m_geom_model`'s entities of
         * dimension `>= d` (never lower, which would over-constrain — e.g. an edge is never pinned
         * to a single point), stopping at the lowest such dimension that has at least one entity
         * within tolerance; every entity within tolerance at that dimension (not just the nearest)
         * is collected into the cell's `geom_targets` — a node near a corner shared by 2+ curves
         * ends up with 2+ targets, by design (see `CellData.h`). The single nearest of those targets
         * is used to snap/refit the cell's own stored geometry: a node's position is projected onto
         * it; a curved edge's interior control points (its 2 endpoints stay pinned to their own,
         * already-classified nodes) are projected onto it; a face's/block's stored surface/volume is
         * then rebuilt from its (now possibly curved) boundary edges/faces via
         * `coons_surface_from_edges()`/`tfi_volume_from_faces()` — never left as a stale blend of the
         * pre-classification straight geometry. Not incremental: every cell is reclassified and
         * every edge/face/block geometry rebuilt from scratch on every call.
         *
         * @param ATolVertex Tolerance for snapping onto a vertex (dimension-0 entity) — the tight
         *        threshold, since 2 distinct vertices are typically much closer to each other than
         *        to the nearest curve/surface.
         * @param ATolCurveSurface Tolerance for snapping onto a curve/surface/volume (dimension >=
         *        1) — typically looser than `ATolVertex`. Defaults to `ATolVertex` when left at its
         *        default `-1.0`.
         */
        void classify(double ATolVertex, double ATolCurveSurface = -1.0) {
            const double tol_cs = (ATolCurveSurface < 0.0) ? ATolVertex : ATolCurveSurface;

            for (auto it = m_cmap.template attributes<0>().begin(), itend = m_cmap.template attributes<0>().end();
                 it != itend;
                 ++it) {
                const auto result = classify_position(GroupDim::Dim0, it->info().point, ATolVertex, tol_cs);
                it->info().geom_targets = result.targets;
                if (result.any()) {
                    project_onto(result.nearest_dim, result.nearest_tag, it->info().point);
                }
            }

            for (auto it = m_cmap.template attributes<1>().begin(), itend = m_cmap.template attributes<1>().end();
                 it != itend;
                 ++it) {
                refit_edge(it, tol_cs);
            }

            for (auto it = m_cmap.template attributes<2>().begin(), itend = m_cmap.template attributes<2>().end();
                 it != itend;
                 ++it) {
                classify_and_rebuild_face(it, tol_cs);
            }

            for (auto it = m_cmap.template attributes<3>().begin(), itend = m_cmap.template attributes<3>().end();
                 it != itend;
                 ++it) {
                classify_and_rebuild_block(it, tol_cs);
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
         * @param ASubdivisions Number of intervals to subdivide every parametric axis into; must be
         *        >= 1. Uniform across the whole blocking (a documented V1 limitation — no
         *        per-block/per-edge subdivision counts, no non-conformal adjacency).
         * @return The generated mesh.
         */
        UnstructuredMesh<CubicTraits> to_mesh(SizeT ASubdivisions) {
            assert(ASubdivisions >= 1 && "Blocking::to_mesh: ASubdivisions must be >= 1");
            const std::size_t s = ASubdivisions;
            UnstructuredMesh<CubicTraits> mesh;

            std::map<Node, NodeId> node_ids;
            for (auto it = m_cmap.template attributes<0>().begin(), itend = m_cmap.template attributes<0>().end();
                 it != itend;
                 ++it) {
                node_ids[it] = mesh.add_node(it->info().point);
            }

            std::map<Edge, EdgeChain> edge_chains;
            for (auto it = m_cmap.template attributes<1>().begin(), itend = m_cmap.template attributes<1>().end();
                 it != itend;
                 ++it) {
                edge_chains[it] = build_edge_chain(it, node_ids, mesh, s);
            }

            std::map<Face, FaceGrid> face_grids;
            for (auto it = m_cmap.template attributes<2>().begin(), itend = m_cmap.template attributes<2>().end();
                 it != itend;
                 ++it) {
                FaceGrid fg = build_face_grid(it, node_ids, edge_chains, mesh, s);
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
                emit_hex_cells(it, node_ids, edge_chains, face_grids, mesh, s);
            }

            return mesh;
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
                    if (positions_match_reversed(dA, dB, 2)) {
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
         * @param ATolVertex Tolerance used only when scanning vertices (dimension 0).
         * @param ATolCurveSurface Tolerance used when scanning curves/surfaces/volumes (dimension >= 1).
         * @return The classification outcome (possibly empty, if nothing is within tolerance anywhere).
         */
        ClassifyResult
        classify_position(GroupDim AMinDim, const Point3d &AP, double ATolVertex, double ATolCurveSurface) const {
            ClassifyResult result;
            double best = 0.0;
            if (AMinDim <= GroupDim::Dim0) {
                scan_dimension(m_geom_model->vertices(), GroupDim::Dim0, AP, ATolVertex, result, best);
                if (result.any()) return result;
            }
            if (AMinDim <= GroupDim::Dim1) {
                scan_dimension(m_geom_model->curves(), GroupDim::Dim1, AP, ATolCurveSurface, result, best);
                if (result.any()) return result;
            }
            if (AMinDim <= GroupDim::Dim2) {
                scan_dimension(m_geom_model->surfaces(), GroupDim::Dim2, AP, ATolCurveSurface, result, best);
                if (result.any()) return result;
            }
            if (AMinDim <= GroupDim::Dim3) {
                scan_dimension(m_geom_model->volumes(), GroupDim::Dim3, AP, ATolCurveSurface, result, best);
            }
            return result;
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
         * @brief Classifies one edge against `m_geom_model` (dimension >= 1 only) and refits its
         * curve: its 2 endpoints are pinned to their (already-classified) nodes' current positions,
         * and its interior control points are re-derived by straight-line interpolation between
         * those endpoints, then projected onto the nearest classification target, if any.
         * @param AEdge The edge to classify and refit.
         * @param ATolCurveSurface Tolerance for the (dimension >= 1) search.
         */
        void refit_edge(Edge AEdge, double ATolCurveSurface) {
            const Dart d = AEdge->dart();
            const Node n0 = m_cmap.template attribute<0>(d);
            const Node n1 = m_cmap.template attribute<0>(m_cmap.template beta<1>(d));
            const Point3d &p0 = n0->info().point;
            const Point3d &p1 = n1->info().point;

            const auto result =
                classify_position(GroupDim::Dim1, AEdge->info().curve.value(0.5), ATolCurveSurface, ATolCurveSurface);
            AEdge->info().geom_targets = result.targets;

            constexpr std::size_t n = TEdgeCurve::NumControlPoints;
            const Vector3d chord(p0, p1);
            TEdgeCurve fitted;
            for (std::size_t i = 0; i < n; ++i) {
                const double t = (n > 1) ? static_cast<double>(i) / static_cast<double>(n - 1) : 0.0;
                Point3d cp = p0 + chord * t;
                if (result.any() && i > 0 && i + 1 < n) {
                    project_onto(result.nearest_dim, result.nearest_tag, cp);
                }
                fitted[i] = cp;
            }
            AEdge->info().curve = fitted;
        }

        /**
         * @brief Builds a control-point-reversed copy of a curve.
         * @param ACurve The curve to reverse.
         * @return The reversed curve.
         */
        static TEdgeCurve reversed_curve(const TEdgeCurve &ACurve) {
            TEdgeCurve rev;
            constexpr std::size_t n = TEdgeCurve::NumControlPoints;
            for (std::size_t i = 0; i < n; ++i) {
                rev[i] = ACurve[n - 1 - i];
            }
            return rev;
        }

        /**
         * @brief Returns @p ACurve as-is if it already starts at @p AExpectedStart, or its reversed
         * copy otherwise — used to orient an already-fitted edge curve consistently for
         * `coons_surface_from_edges()`/`tfi_volume_from_faces()` regardless of which of its 2
         * incident faces/blocks last stored it.
         * @param ACurve The curve to orient.
         * @param AExpectedStart The position its first control point must have.
         * @return @p ACurve, oriented to start at @p AExpectedStart.
         */
        static TEdgeCurve oriented_curve(const TEdgeCurve &ACurve, const Point3d &AExpectedStart) {
            if (ACurve.control_points()[0] == AExpectedStart) return ACurve;
            return reversed_curve(ACurve);
        }

        /**
         * @brief Classifies one face against `m_geom_model` (dimension >= 2 only, using the face's
         * own corner centroid) and rebuilds its stored surface via `coons_surface_from_edges()` from
         * its 4 (now possibly refitted) boundary edges. Applies uniformly to a standalone 2D (quad)
         * block's own face and to one of a 3D (hex) block's 6 bounding faces — both are ordinary
         * 2-cells with exactly 4 boundary edges.
         * @param AFace The face to classify and rebuild.
         * @param ATolCurveSurface Tolerance for the (dimension >= 2) search.
         */
        void classify_and_rebuild_face(Face AFace, double ATolCurveSurface) {
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
            const auto result = classify_position(GroupDim::Dim2, center, ATolCurveSurface, ATolCurveSurface);
            AFace->info().geom_targets = result.targets;

            std::array<TEdgeCurve, 4> curves{};
            for (auto it = m_cmap.template one_dart_per_incident_cell<1, 2>(fd).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<1, 2>(fd).end();
                 it != itend;
                 ++it) {
                const int a = node_index(local_nodes, m_cmap.template attribute<0>(it));
                const int b = node_index(local_nodes, m_cmap.template attribute<0>(m_cmap.template beta<1>(it)));
                const std::size_t k = find_edge(QUAD_EDGES, a, b);
                const Edge e = m_cmap.template attribute<1>(it);
                curves[k] = oriented_curve(e->info().curve,
                                           local_nodes[static_cast<std::size_t>(QUAD_EDGES[k].first)]->info().point);
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
         * @param ABlock The block to classify and rebuild.
         * @param ATolCurveSurface Tolerance for the (dimension >= 3) search.
         */
        void classify_and_rebuild_block(Block ABlock, double ATolCurveSurface) {
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
            const auto result = classify_position(GroupDim::Dim3, center, ATolCurveSurface, ATolCurveSurface);
            ABlock->info().geom_targets = result.targets;

            std::array<TEdgeCurve, 12> curves{};
            for (auto it = m_cmap.template one_dart_per_incident_cell<1, 3>(bd).begin(),
                      itend = m_cmap.template one_dart_per_incident_cell<1, 3>(bd).end();
                 it != itend;
                 ++it) {
                const int a = node_index(local_nodes, m_cmap.template attribute<0>(it));
                const int b = node_index(local_nodes, m_cmap.template attribute<0>(m_cmap.template beta<1>(it)));
                const std::size_t k = find_edge(HEX_EDGES, a, b);
                const Edge e = m_cmap.template attribute<1>(it);
                curves[k] = oriented_curve(e->info().curve,
                                           local_nodes[static_cast<std::size_t>(HEX_EDGES[k].first)]->info().point);
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
         * @brief Builds one edge's `to_mesh()` node-id chain: its 2 endpoints reuse the already
         * mapped corner node ids, its `S-1` interior points are newly added mesh nodes sampled
         * directly from the edge's own curve.
         * @param AEdge The edge to sample.
         * @param ANodeIds Already-populated node-attribute-to-mesh-id map (its corner endpoints).
         * @param AMesh The mesh being built, appended to for interior points.
         * @param AS Subdivisions per axis.
         * @return The built chain.
         */
        EdgeChain build_edge_chain(Edge AEdge,
                                   const std::map<Node, NodeId> &ANodeIds,
                                   UnstructuredMesh<CubicTraits> &AMesh,
                                   std::size_t AS) {
            const Dart d = AEdge->dart();
            const Node na = m_cmap.template attribute<0>(d);
            const Node nb = m_cmap.template attribute<0>(m_cmap.template beta<1>(d));
            const bool start_is_a = (AEdge->info().curve.control_points()[0] == na->info().point);

            EdgeChain chain;
            chain.start = start_is_a ? na : nb;
            const Node end = start_is_a ? nb : na;
            chain.ids.resize(AS + 1);
            chain.ids[0] = ANodeIds.at(chain.start);
            chain.ids[AS] = ANodeIds.at(end);
            for (std::size_t i = 1; i < AS; ++i) {
                const double t = static_cast<double>(i) / static_cast<double>(AS);
                chain.ids[i] = AMesh.add_node(AEdge->info().curve.value(t));
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
         * @param AS Subdivisions per axis.
         * @return The built grid.
         */
        FaceGrid build_face_grid(Face AFace,
                                 const std::map<Node, NodeId> &ANodeIds,
                                 const std::map<Edge, EdgeChain> &AEdgeChains,
                                 UnstructuredMesh<CubicTraits> &AMesh,
                                 std::size_t AS) {
            FaceGrid fg;
            const Dart fd = AFace->dart();
            Dart walk = fd;
            for (std::size_t c = 0; c < 4; ++c) {
                fg.local_nodes[c] = m_cmap.template attribute<0>(walk);
                walk = m_cmap.template beta<1>(walk);
            }

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
         * @param AS Subdivisions per axis.
         */
        void emit_hex_cells(Block ABlock,
                            const std::map<Node, NodeId> &ANodeIds,
                            const std::map<Edge, EdgeChain> &AEdgeChains,
                            const std::map<Face, FaceGrid> &AFaceGrids,
                            UnstructuredMesh<CubicTraits> &AMesh,
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
        static TEdgeCurve straight_curve(const Point3d &AA, const Point3d &AB) {
            TEdgeCurve curve;
            constexpr std::size_t n = TEdgeCurve::NumControlPoints;
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

        /** @brief Non-owning pointer to the geometric model this blocking is built against. */
        const TGeomModel *m_geom_model;
        /** @brief The underlying (always dimension-3) combinatorial map. */
        Map m_cmap;
        /** @brief Every face created as one of a hex block's 6 boundary faces (dim-3 sew candidates). */
        std::vector<Face> m_hex_faces;
    };

} // namespace gecko
