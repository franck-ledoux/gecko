#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include <gecko/block/Blocking.h>
#include <gecko/math/BezierCurve.h>
#include <gecko/math/Point3d.h>

#include <gecko/app/GeomModelFacade.h>

namespace gecko::app {

    /**
     * @class BlockingFacade
     * @brief Python-facing façade over gecko::Blocking<FacetedGeometry, TEdgeCurve>: every method
     * takes/returns only primitive types — corners as lists of (x,y,z) tuples, faces/blocks as
     * plain int ids assigned by this façade, never a gecko::Point3d or CGAL dart/attribute handle
     * — see docs/user-guide/python.md.
     *
     * The edge curve degree (1 = straight, higher = Bezier of that order) is a plain int, given at
     * construction and changeable afterwards through set_degree(). It used to be neither: the degree
     * was a C++ template parameter, so this façade held a std::variant of one implementation per
     * degree, could offer only the degrees it had instantiated, and could not change one without
     * building a whole new Blocking.
     */
    class BlockingFacade {
    public:
        /** @brief Lowest usable edge curve degree. 1 is a straight edge. */
        static constexpr int MIN_DEGREE = 1;

        /**
         * @brief Constructor.
         * @param model The geometric model to build against; must outlive this Blocking (enforced
         *        on the Python side via py::keep_alive).
         * @param degree Edge curve degree: 1 for straight edges, higher for Bezier curves of that
         *        order. Changeable afterwards through set_degree().
         * @throw std::invalid_argument if @p degree is below MIN_DEGREE.
         */
        explicit BlockingFacade(const GeomModelFacade &model, int degree = 1);

        /** @brief Gets the edge curve degree this blocking currently uses. @return The degree. */
        [[nodiscard]] int degree() const;

        /**
         * @brief Rebuilds every edge, face and block at a new degree, refitting them onto the model.
         *
         * Topology and classification are untouched; only the representation changes. Raising the
         * order does not merely add control points, it uses them: an edge lying on a curved model
         * curve, which at degree 1 could only be its chord, comes back following it.
         * @param degree The new degree, at least MIN_DEGREE.
         * @param tol_vertex Tolerance for snapping onto a vertex, as classify() takes it.
         * @param tol_curve Tolerance for snapping onto a curve; defaults to @p tol_vertex.
         * @param tol_surface Tolerance for snapping onto a surface; defaults to the resolved curve
         *        tolerance.
         * @throw std::invalid_argument if @p degree is below MIN_DEGREE.
         */
        void set_degree(int degree, double tol_vertex, double tol_curve = -1.0, double tol_surface = -1.0);

        /**
         * @brief Creates a new, unsewn quad block.
         * @param corners The 4 corner positions, each an (x,y,z) triple, in perimeter order.
         * @return The id assigned to the newly created face.
         * @throw std::invalid_argument if @p corners doesn't have exactly 4 entries.
         */
        int create_quad_block(const std::vector<std::array<double, 3>> &corners);

        /**
         * @brief Creates a new, unsewn hex block.
         * @param corners The 8 corner positions, each an (x,y,z) triple — see
         *        gecko::Blocking::create_hex_block for the expected HEX8 ordering.
         * @return The id assigned to the newly created block, which keeps naming it for as long as
         *         it exists.
         * @throw std::invalid_argument if @p corners doesn't have exactly 8 entries.
         */
        int create_hex_block(const std::vector<std::array<double, 3>> &corners);

        /**
         * @brief How far each edge departs from a straight line: the largest distance from one of its
         * interior control points to the chord joining its 2 endpoints, per edge.
         *
         * A diagnostic rather than a modelling quantity. A straight blocking must read as zero here
         * whatever has been done to it — cutting a straight edge gives straight halves, and deleting
         * moves nothing — so a non-zero entry with nothing classified to justify it means the
         * geometry has been corrupted, and a zero one while the screen shows curves means the
         * drawing has.
         *
         * @return One deviation per edge, in the order edge_vertices()/edge_classification_dims()
         *         use. Always 0 at degree 1, whose edges cannot bend.
         */
        [[nodiscard]] std::vector<double> edge_bends() const;

        /**
         * @brief Deletes one block, along with every face, edge and corner that existed only because
         * of it — what it shared with a neighbouring block stays, as that neighbour's boundary.
         * @param block_id A block id, as block_ids() reports them.
         * @throw std::out_of_range if @p block_id is not a block of this blocking.
         */
        void delete_block(int block_id);

        /** @brief Auto-detects and sews coincident blocks created so far. Not incremental. */
        void build_connectivity();

        /**
         * @brief Classifies every cell onto the geometric model and refits geometry accordingly.
         * @param tol_vertex Tolerance for snapping onto a vertex.
         * @param tol_curve Tolerance for snapping onto a curve; defaults to @p tol_vertex.
         * @param tol_surface Tolerance for snapping onto a surface (and volume); defaults to the
         *        resolved curve tolerance.
         */
        void classify(double tol_vertex, double tol_curve = -1.0, double tol_surface = -1.0);

        /**
         * @brief Snaps one corner node onto the geometric model, reclassifying and refitting only
         * the cells that touch it — meant to run when a dragged corner is released.
         * @param node_id A node id from node_ids().
         * @param tol_vertex Tolerance for snapping onto a vertex.
         * @param tol_curve Tolerance for snapping onto a curve; defaults to @p tol_vertex.
         * @param tol_surface Tolerance for snapping onto a surface; defaults to the resolved curve
         *        tolerance.
         * @throw std::out_of_range if @p node_id is not a known node id.
         */
        void snap_node(int node_id, double tol_vertex, double tol_curve = -1.0, double tol_surface = -1.0);

        /**
         * @brief Gets the number of cells of a given dimension in the blocking.
         * @param dim Cell dimension, in [0,3].
         * @return The cell count.
         * @throw std::invalid_argument if @p dim is not in [0,3].
         */
        [[nodiscard]] std::size_t nb_cells(int dim) const;

        /** @brief Checks the topological validity of the underlying map. @return true if valid. */
        [[nodiscard]] bool is_valid_topology() const;
        /** @brief Checks whether this blocking has no 3-cell (hex block). @return true if purely 2D. */
        [[nodiscard]] bool is_purely_2d() const;
        /**
         * @brief Checks whether a face can be deleted.
         * @param face_id A face id returned by create_quad_block() (or a hex block's boundary
         *        face — not currently reachable through this façade).
         * @return true if the face can be deleted.
         * @throw std::out_of_range if @p face_id is not a known face id.
         */
        [[nodiscard]] bool can_delete_face(int face_id) const;
        /**
         * @brief Deletes a face from the structure.
         * @param face_id A face id returned by create_quad_block().
         * @throw std::out_of_range if @p face_id is not a known face id.
         */
        void delete_face(int face_id);

        /**
         * @brief Generates a quad/hex mesh reproducing the blocking's geometry and writes it to a
         * VTK legacy ASCII file.
         * @param subdivisions Number of intervals to subdivide every parametric axis into (>= 1).
         * @param path Output file path.
         */
        void write_vtk(int subdivisions, const std::string &path);

        /** @brief Gets the id of every corner node of the block structure. @return The node ids. */
        [[nodiscard]] std::vector<int> node_ids() const;

        /**
         * @brief Every edge's id, in the order the display accessors of that dimension emit —
         * edge_vertices(), edge_segments(), edge_bends(), edge_classification_dims().
         *
         * The bridge between the two ways of naming a cell here, and both are needed. What is drawn
         * has to be positional: a renderer indexes a flat array. What is *acted on* has to be an id:
         * a position is only true until the next operation renumbers it, and a viewer that picks a
         * position one frame and acts on it the next would act on a different edge.
         *
         * @return One id per edge, aligned with those accessors.
         */
        [[nodiscard]] std::vector<int> edge_ids() const;
        /** @brief Every face's id, in the order face_classification_dims() and face_grid_owners()
         * index. @return One id per face. @see edge_ids() */
        [[nodiscard]] std::vector<int> face_ids() const;
        /** @brief Every block's id, in the order block_volumes() and mesh_hex_owners() index.
         * @return One id per block. @see edge_ids() */
        [[nodiscard]] std::vector<int> block_ids() const;

        /**
         * @brief The 6 faces bounding one block.
         *
         * What a nappe is named with: pillow() takes a set of face ids, and the faces of a block are
         * where any nappe closed around something starts from.
         *
         * @param block_id A block id, as block_ids() reports them.
         * @return Its 6 face ids, in no particular order.
         * @throw std::out_of_range if no block carries that id.
         */
        [[nodiscard]] std::vector<int> block_faces(int block_id);

        /**
         * @brief The blocks one face bounds — 2 of them, or 1 where the face is on the boundary of
         * the blocking.
         *
         * The other half of naming a nappe: which side of a face is which, so that a caller can say
         * which of the 2 pillow() should shrink.
         *
         * @param face_id A face id, as face_ids() reports them.
         * @return The ids of the blocks it bounds — empty for a standalone quad block, which bounds
         *         none.
         * @throw std::out_of_range if no face carries that id.
         */
        [[nodiscard]] std::vector<int> face_blocks(int face_id);
        /**
         * @brief Gets the position of a corner node.
         * @param node_id A node id, from node_ids().
         * @return Its (x,y,z) position.
         * @throw std::out_of_range if @p node_id is not a known node id.
         */
        [[nodiscard]] std::array<double, 3> node_position(int node_id) const;
        /**
         * @brief Moves a corner node, refitting every incident edge/face/block geometry.
         * @param node_id A node id, from node_ids().
         * @param x New X coordinate.
         * @param y New Y coordinate.
         * @param z New Z coordinate.
         * @throw std::out_of_range if @p node_id is not a known node id.
         */
        void move_node(int node_id, double x, double y, double z);

        /**
         * @brief Gets what each corner node is classified on, in node_ids() order.
         * @return One entry per node: the dimension of the geometric entity it is classified on
         *         (0 vertex, 1 curve, 2 surface, 3 volume), or -1 if it is unclassified. Every
         *         node is unclassified until classify() runs.
         */
        [[nodiscard]] std::vector<int> node_classification_dims() const;

        /**
         * @brief Gets what each block edge is classified on, in the traversal order edge_vertices()
         * and edge_segments() use.
         * @return One entry per edge: the dimension of the geometric entity it is classified on
         *         (0 vertex, 1 curve, 2 surface, 3 volume), or -1 if unclassified.
         */
        [[nodiscard]] std::vector<int> edge_classification_dims() const;

        /**
         * @brief Gets what each block face is classified on, in the traversal order
         * face_grid_owners() indexes into.
         * @return One entry per face: the dimension of the geometric entity it is classified on, or
         *         -1 if unclassified.
         */
        [[nodiscard]] std::vector<int> face_classification_dims() const;

        /**
         * @brief Gets every block edge's control points — the handles actually driving a curved
         * edge's shape, as opposed to the points sampled along it by edge_vertices().
         *
         * Only informative above degree 1, where a straight edge's 2 control points are just its
         * endpoints.
         * @return `degree + 1` (x,y,z) triples per edge, edge after edge, in the same traversal
         *         order as edge_vertices().
         */
        [[nodiscard]] std::vector<std::array<double, 3>> edge_control_points() const;

        /**
         * @brief The segments of every edge's control polygon, joining consecutive control points
         * within each edge (never across two edges).
         * @return One pair of edge_control_points() indices per segment.
         */
        [[nodiscard]] std::vector<std::array<int, 2>> edge_control_polygons() const;

        /**
         * @brief Gets every block face's control points — the `(degree+1)²` grid driving its Bezier
         * surface, the 2D counterpart of edge_control_points().
         * @return `(degree+1)²` (x,y,z) triples per face, face after face, row-major (`u` outer,
         *         `v` inner), in the traversal order face_classification_dims() uses.
         */
        [[nodiscard]] std::vector<std::array<double, 3>> face_control_points() const;
        /**
         * @brief The segments of every face's control net, joining each control point to its
         * neighbours along `u` and `v` within one face.
         * @return One pair of face_control_points() indices per segment.
         */
        [[nodiscard]] std::vector<std::array<int, 2>> face_control_nets() const;

        /**
         * @brief Gets every block's control points — the `(degree+1)³` grid driving its Bezier
         * volume, the 3D counterpart of edge_control_points().
         * @return `(degree+1)³` (x,y,z) triples per block, block after block, row-major
         *         (`u` outermost, `w` innermost).
         */
        [[nodiscard]] std::vector<std::array<double, 3>> block_control_points() const;
        /**
         * @brief The segments of every block's control lattice, joining each control point to its
         * neighbours along `u`, `v` and `w` within one block.
         * @return One pair of block_control_points() indices per segment.
         */
        [[nodiscard]] std::vector<std::array<int, 2>> block_control_lattices() const;

        /**
         * @brief Samples every block face into a grid of quads, for display.
         *
         * Unlike mesh_quads(), which only emits anything for standalone 2D blocks, this covers the
         * bounding faces of 3D blocks too — so a blocking's faces can be drawn (and colored by
         * classification) whatever its dimension. Use with face_grid_quads() and face_grid_owners().
         * @param subdivisions Number of intervals per parametric axis (>= 1).
         * @return `(subdivisions+1)^2` (x,y,z) triples per face, face after face.
         */
        [[nodiscard]] std::vector<std::array<double, 3>> face_grid_vertices(int subdivisions) const;
        /**
         * @brief The quads joining the points face_grid_vertices() returns.
         * @param subdivisions Number of intervals per parametric axis (>= 1).
         * @return One quadruple of face_grid_vertices() indices per quad.
         */
        [[nodiscard]] std::vector<std::array<int, 4>> face_grid_quads(int subdivisions) const;
        /**
         * @brief Which block face each quad of face_grid_quads() came from, so a per-face value
         * (its classification, say) can be spread onto every quad subdividing it.
         * @param subdivisions Number of intervals per parametric axis (>= 1).
         * @return One index into face_classification_dims() per quad.
         */
        [[nodiscard]] std::vector<int> face_grid_owners(int subdivisions) const;

        /**
         * @brief Samples every edge of the block structure along its own curve.
         *
         * Together with edge_segments(), gives the block structure's own edges as polylines —
         * distinct from the subdivision lines of the mesh the blocking generates. For a curved
         * (degree 3) blocking, @p samples > 1 traces the actual curve rather than its chord.
         * @param samples Number of intervals to split each edge into (>= 1).
         * @return One (x,y,z) triple per sample point, `samples + 1` per edge, edge after edge.
         */
        [[nodiscard]] std::vector<std::array<double, 3>> edge_vertices(int samples) const;
        /**
         * @brief The segments joining the points edge_vertices() returns, in the same traversal
         * order (call the two together, without mutating the blocking in between).
         * @param samples Number of intervals to split each edge into (>= 1).
         * @return One pair of edge_vertices() indices per segment.
         */
        [[nodiscard]] std::vector<std::array<int, 2>> edge_segments(int samples) const;

        /**
         * @brief Generates the mesh reproducing the blocking and returns its node positions.
         * @param subdivisions Number of intervals to subdivide every parametric axis into (>= 1).
         * @return One (x,y,z) triple per mesh node.
         */
        [[nodiscard]] std::vector<std::array<double, 3>> mesh_vertices(int subdivisions);
        /**
         * @brief Generates the mesh reproducing the blocking and returns its quads.
         * @param subdivisions Number of intervals to subdivide every parametric axis into (>= 1).
         * @return One quadruple of mesh_vertices() indices per quad.
         */
        [[nodiscard]] std::vector<std::array<int, 4>> mesh_quads(int subdivisions);
        /**
         * @brief Generates the mesh reproducing the blocking and returns its hexes.
         * @param subdivisions Number of intervals to subdivide every parametric axis into (>= 1).
         * @return One 8-tuple of mesh_vertices() indices per hex.
         */
        [[nodiscard]] std::vector<std::array<int, 8>> mesh_hexes(int subdivisions);
        /**
         * @brief Which block each hex of mesh_hexes() was generated from, so a per-block value can
         * be spread onto every cell subdividing it — one colour per block, say.
         * @param subdivisions Number of intervals per parametric axis (>= 1); must match the one
         *        mesh_hexes() was called with, since it sets how many cells each block contributes.
         * @return One index into the block order — the one block_volumes()/delete_block() speak —
         *         per entry of mesh_hexes().
         */
        [[nodiscard]] std::vector<int> mesh_hex_owners(int subdivisions);
        /**
         * @brief Which block each quad of mesh_quads() was generated from — the 2D counterpart of
         * mesh_hex_owners().
         * @param subdivisions Number of intervals per parametric axis (>= 1); must match the one
         *        mesh_quads() was called with.
         * @return One index per entry of mesh_quads(), counting standalone quad blocks in the order
         *         they are emitted in. A hex's own bounding faces generate no mesh quads and are not
         *         counted, so this is *not* an index into face_classification_dims().
         */
        [[nodiscard]] std::vector<int> mesh_quad_owners(int subdivisions);

        /**
         * @brief Measures every block of the blocking, from its own stored geometry.
         * @param subdivisions Intervals per parametric axis (>= 1). Exact at 1 for a block whose
         *        faces are planar; converges as it grows for a curved or warped one.
         * @return One signed volume per block, in the order block_ids() reports them. A negative
         *         value means that block's frame is inverted.
         */
        [[nodiscard]] std::vector<double> block_volumes(int subdivisions);

        /**
         * @brief The edges cut_sheet() would split if aimed at @p edge_id — the whole sheet, not just
         * the one edge — so a caller can show what a cut is about to do before doing it.
         * @param edge_id An edge id, as edge_ids() reports them.
         * @return A *position* per edge of the sheet, in the order edge_vertices()/edge_segments()
         *         use, or an empty list when the sheet cannot be cut homogeneously. Positions rather
         *         than ids because this answer exists to be drawn: see edge_ids().
         * @throw std::out_of_range if @p edge_id is not an edge of this blocking.
         */
        [[nodiscard]] std::vector<int> sheet_edges(int edge_id);

        /**
         * @brief Where a cut would actually land: one point per sheet edge, each on the side the
         * sheet agrees on, so the whole cut locus can be shown before committing to it.
         * @param edge_id An edge, in the same order sheet_edges() uses.
         * @param param Where along that edge to cut, strictly inside (0, 1).
         * @return One point per edge of the sheet, in sheet_edges() order, or an empty list when the
         *         sheet cannot be cut homogeneously.
         * @throw std::out_of_range if @p edge_id is not an edge of this blocking.
         */
        [[nodiscard]] std::vector<std::array<double, 3>> sheet_cut_points(int edge_id, double param);

        /**
         * @brief Cuts the blocking along the whole sheet through one edge, splitting every block the
         * sheet crosses in 2 and keeping the geometry it cut through exactly (see
         * `Blocking::cut_sheet()`).
         * @param edge_id An edge, in the same order sheet_edges() uses.
         * @param param Where along that edge to cut, measured along its own curve exactly as
         *        edge_vertices() samples it, strictly inside (0, 1).
         * @return false, changing nothing, if @p param is out of range or the sheet cannot be cut
         *         homogeneously.
         * @throw std::out_of_range if @p edge_id is not an edge of this blocking.
         */
        bool cut_sheet(int edge_id, double param);

        /**
         * @brief Deletes the whole sheet through one edge, gluing back what was either side of it —
         * the inverse of cut_sheet() (see `Blocking::delete_sheet()`).
         *
         * Merging 2 corners into 1 means one of their 2 classifications has to go, and the more
         * constrained side wins: a corner on a model vertex beats one on a curve, which beats one on
         * a surface, and the merged corner is projected onto the winner. So a block structure fitted
         * to a model does not drift off its features when a layer is taken out of it.
         *
         * @param edge_id An edge of the sheet, in the same order sheet_edges() uses.
         * @param tol_vertex Tolerance for snapping onto a vertex, as classify() defines it — used
         *        only where a refitted cell has to fall back on a proximity search.
         * @param tol_curve Tolerance for snapping onto a curve. Defaults to @p tol_vertex.
         * @param tol_surface Tolerance for snapping onto a surface. Defaults to the curve one.
         * A sheet holding every block there is collapses like any other and leaves the blocking
         * empty — a state to be in rather than a broken one, and what taking the last layer out
         * means. The count of blocks is never a reason to refuse; the geometry is the only thing that
         * ever stands in the way, and it does so in exactly one case: an edge of the sheet joining 2
         * corners classified on 2 *different* model vertices. Merging those would leave one of the 2
         * vertices with no corner of the block structure on it, and nothing else records where it
         * was.
         *
         * @return false, changing nothing, when the sheet cannot be collapsed — see
         *         `Blocking::delete_sheet()` for the cases.
         * @throw std::out_of_range if @p edge_id is not an edge of this blocking.
         */
        bool delete_sheet(int edge_id, double tol_vertex, double tol_curve = -1.0, double tol_surface = -1.0);

        /**
         * @brief Inserts a layer of blocks along a nappe of block faces — the pillowing operation
         * (see `Blocking::pillow()`).
         *
         * The nappe is a sheet of faces that cuts the blocking in two: closed around a set of blocks
         * to isolate them, or running clean through the structure and out on its boundary. Which of
         * the 2 sides shrinks has to be said, @p inside_block_id saying it — the layer is inserted
         * into the gap that side leaves, and the other side does not move at all. Where the nappe
         * lies on the boundary of the blocking, that other side is the model's own boundary, and
         * this is what keeps the structure on the geometry it was classified onto.
         *
         * A corner the nappe cuts through becomes 2. The outside one keeps its classification; the
         * inside one keeps only what it is still on after moving, and comes back unclassified when
         * it has been pushed off everything. Nothing is classified by proximity here, so a blocking
         * nobody classified stays that way.
         *
         * @param face_ids The nappe, each face named once.
         * @param inside_block_id A block on the side that shrinks.
         * @param thickness How far that side is pulled back, as a fraction of the mean length of the
         *        edges at each corner that moves. In `(0,1)`.
         * @param tol_vertex Tolerance for a moved corner staying on a vertex, as classify() defines it.
         * @param tol_curve Tolerance for staying on a curve. Defaults to @p tol_vertex.
         * @param tol_surface Tolerance for staying on a surface. Defaults to the curve one.
         * @return false, changing nothing, when what was given is not a nappe — see
         *         `Blocking::pillow()` for the cases.
         * @throw std::out_of_range if one of @p face_ids is not a face of this blocking, or
         *        @p inside_block_id not a block of it.
         */
        bool pillow(const std::vector<int> &face_ids,
                    int inside_block_id,
                    double thickness,
                    double tol_vertex,
                    double tol_curve = -1.0,
                    double tol_surface = -1.0);

        /**
         * @brief Whether there is an edit to take back. @return true if undo() would do something.
         */
        [[nodiscard]] bool can_undo() const;
        /** @brief Whether there is an undone edit to put back. @return true if redo() would do
         * something. */
        [[nodiscard]] bool can_redo() const;

        /**
         * @brief Takes back the most recent edit, restoring the blocking as it was just before it.
         *
         * Snapshot-based: every operation that changes the blocking copies it first, and undoing puts
         * that copy back. Not a stack of inverse operations, because the operations here do not have
         * inverses — collapsing the layer a cut just created does not restore the block, as
         * `delete_sheet()` documents and the tests fix.
         *
         * Ids survive an undo, since the snapshot carries them: a caller holding a block id from
         * before the undone edit finds its block again, exactly as it would had the edit not happened.
         *
         * Does nothing when there is nothing to take back.
         */
        void undo();
        /** @brief Puts back the most recently undone edit. Does nothing when there is none; any new
         * edit discards the redo history, which is what makes the two stacks a line rather than a
         * tree. */
        void redo();

        /**
         * @brief How many edits can be taken back. Beyond this the oldest is dropped.
         *
         * A snapshot is a whole copy of the blocking, so the history costs depth times its size.
         * A few thousand cells at degree 3 is a couple of megabytes each, which is what makes 20 a
         * sensible default rather than a compromise.
         *
         * @param depth Number of edits to keep, at least 1.
         * @throw std::invalid_argument if @p depth is below 1.
         */
        void set_history_depth(int depth);
        /** @brief How many edits can be taken back. @return The current depth. */
        [[nodiscard]] int history_depth() const;

        /** @brief The blocking and what little bookkeeping is still kept alongside it; public only
         * so the .cpp's free-function helpers can name it — never part of the class' actual
         * (Python-facing) interface. */
        struct Impl {
            /** @brief The kernel type this façade wraps. */
            using BlockingT = Blocking<FacetedGeometry>;

            /** @brief The blocking itself — all the state there is, now that cells name themselves. */
            BlockingT blocking;

            /**
             * @brief Builds the blocking this façade wraps.
             * @param geom The geometric model to build against.
             * @param degree The Bezier degree every cell's geometry is built at.
             */
            explicit Impl(const FacetedGeometry &geom, int degree) : blocking(geom, static_cast<std::size_t>(degree)) {}
        };

    private:
        /**
         * @brief Scope guard taking the snapshot one edit can be undone from.
         *
         * Declared at the top of every method that changes the blocking, for the reason
         * `Blocking::EditSession` is: written as a call, "snapshot before you edit" is a rule to
         * remember, and the one place it is forgotten is the one edit that cannot be taken back.
         *
         * An operation that *refuses* — `cut_sheet()` and `delete_sheet()` both can — calls
         * `discard()`, since a snapshot of a state nothing changed would spend an undo step undoing
         * nothing.
         */
        class Checkpoint {
        public:
            /** @brief Snapshots @p AFacade's blocking. @param AFacade The façade being edited. */
            explicit Checkpoint(BlockingFacade &AFacade);
            Checkpoint(const Checkpoint &) = delete;
            Checkpoint &operator=(const Checkpoint &) = delete;
            /** @brief Keeps the snapshot, unless it was discarded. */
            ~Checkpoint();
            /** @brief Drops the snapshot: the edit did not happen. */
            void discard() { m_keep = false; }

        private:
            /** @brief The façade whose history this guard writes into. */
            BlockingFacade &m_facade;
            /** @brief Whether the snapshot is still wanted when the guard ends. */
            bool m_keep = true;
        };

        Impl m_impl;
        /** @brief Snapshots to undo to, oldest first; the last is the state before the last edit. */
        std::vector<Impl::BlockingT> m_undo;
        /** @brief Snapshots to redo to, discarded by any new edit. */
        std::vector<Impl::BlockingT> m_redo;
        /** @brief How many entries `m_undo` keeps. @see set_history_depth() */
        int m_history_depth = 20;
    };

} // namespace gecko::app
