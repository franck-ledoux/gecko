#pragma once

#include <atomic>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

#include <glm/glm.hpp>

namespace polyscope {
    class Structure;
} // namespace polyscope

#include "BiyConfig.h"
#include <gecko/app/BlockingFacade.h>
#include <gecko/app/GeomModelFacade.h>

namespace gecko::biy {

    /**
     * @brief What the left mouse button does. Polyscope processes its own camera navigation
     * *before* the per-frame user callback runs, so a drag can't be intercepted after the fact —
     * navigation has to be switched off ahead of time, which makes this a genuine mode rather than
     * a per-event decision.
     */
    enum class MouseMode {
        Camera,   ///< Polyscope's usual navigation: rotate/pan/zoom the view.
        Edit,     ///< Navigation off; dragging picks up and moves a block corner.
        Cut,      ///< Navigation off; hovering an edge previews a sheet cut, clicking performs it.
        Collapse, ///< Navigation off; hovering an edge lights its sheet, clicking collapses it away.
        Delete,   ///< Navigation off; hovering a block highlights it, clicking deletes it.
        Pillow,   ///< Navigation off; clicking block faces gathers a nappe, Space inserts a layer along it.
        Chord,    ///< Navigation off; hovering a face shows the fold, clicking folds its chord away.
        Open      ///< Navigation off; an edge then 2 faces round it, and a column opens along the chain.
    };

    /**
     * @enum MeshColoring
     * @brief How the generated mesh is coloured.
     *
     * Subdividing a blocking hides the very thing it was built from: past 2 or 3 intervals the block
     * edges disappear under the cells, and a uniformly coloured mesh says nothing about which block
     * any part of it came from. One hue per block puts that back without needing the block structure
     * drawn on top of it.
     */
    enum class MeshColoring {
        Uniform, ///< One `mesh_color` for the whole mesh.
        ByBlock  ///< One hue per block, spread over every cell that block generated.
    };

    /**
     * @class BiyApp
     * @brief Owns biy's whole live state — the geometric model, the blocking being built against
     * it, and the Polyscope structures displaying both — and drives one frame of the UI at a time.
     *
     * The same state is reachable from two threads: this class' own per-frame work on the main
     * (render) thread, and the embedded Python REPL on its background thread (see PythonConsole).
     * `mutex()` guards every access to the facades; the REPL locks it around each statement it
     * executes, and every method here locks it as needed.
     */
    class BiyApp {
    public:
        /** @brief Highest display subdivision the panel accepts. A 3D blocking's generated mesh
         * grows cubically with it — 20 is already 8000 hexes per block — and the window becomes
         * unusable well before anything else complains. */
        static constexpr int MAX_SUBDIVISIONS = 20;

        /** @brief Highest edge order the panel offers. Nothing in the kernel caps the degree — it is
         * a plain number carried by each curve — but a fit has only so much to say: past this the
         * extra control points crowd together without following the model any closer, and each one
         * costs on every face and block that touches the edge. */
        static constexpr int MAX_ORDER = 10;

        /**
         * @brief Constructor. Loads the geometric model and registers it for display.
         * @param model_path Path to the .msh file to load.
         * @param order Edge order of the blocking to build (see `BlockingFacade`'s degree). Fixed
         *        for the session: it is baked into the C++ blocking type.
         * @throw std::runtime_error if the file cannot be read.
         */
        BiyApp(const std::string &model_path, int order);

        /** @brief Draws this frame's ImGui panel and services any in-progress corner drag. Called
         * by Polyscope once per frame, on the main thread. */
        void per_frame();

        /** @brief The lock guarding the facades, shared with the Python console. @return The mutex. */
        [[nodiscard]] std::mutex &mutex() { return m_mutex; }

        /** @brief The geometric model, for injection into the Python console. @return The model. */
        [[nodiscard]] app::GeomModelFacade &model() { return *m_model; }

        /**
         * @brief Gets the blocking, creating an empty one on first use so the Python console always
         * has a `blocking` to talk to.
         * @return The blocking.
         */
        [[nodiscard]] app::BlockingFacade &blocking();

        /**
         * @brief Asks for the Polyscope structures to be rebuilt on the next frame. Call after any
         * mutation made off the render thread — in practice, from the Python console.
         *
         * Deliberately not a direct `refresh_view()`: registering and removing Polyscope structures
         * races with the main thread drawing them. `mutex()` cannot prevent that, since it only
         * guards this class' own per-frame callback and not Polyscope's rendering, which happens
         * outside it. Deferring to the next frame keeps every structure mutation on the render
         * thread, where drawing cannot be interleaved with it.
         *
         * Safe to call with `mutex()` already held — as the console does — since it takes no lock.
         */
        void request_refresh();

    private:
        /** @brief Re-reads the facades and rebuilds every Polyscope structure. Main thread only —
         * see `request_refresh()`. */
        void refresh_view();
        /** @brief Registers the model's own facets (triangles, and tets when the file had any). */
        void register_model();
        /** @brief Freezes the scene's bounding box/length scale to the model, so later edits to
         * the blocking can't move the ground plane or rescale the view. */
        void freeze_scene_extents();
        /** @brief Switches what the left mouse button does, turning Polyscope's own navigation on
         * or off to match. */
        void set_mouse_mode(MouseMode mode);
        /** @brief Shows the drag highlight on one corner, or hides it when @p ANodeId is unset. */
        void show_highlight(std::optional<int> node_id);
        /** @brief Applies a per-cell classification color quantity to a Polyscope structure.
         * @tparam TStructure Polyscope structure type.
         * @tparam TAdd Callable adding a color quantity of the right kind to that structure.
         * @param AStructure The structure to color.
         * @param ADims One classification dimension per element (-1 when unclassified).
         * @param AAdd Adds the quantity (structures differ: per-edge, per-face, ...).
         */
        template<typename TStructure, typename TAdd>
        void apply_classification_colors(TStructure *structure, const std::vector<int> &dims, TAdd add);
        /** @brief Registers (or removes) all three control displays: edge polygons, face nets and
         * block lattices. */
        void refresh_control_nets();
        /**
         * @brief Registers (or removes) one control display: its points, drawn as spheres, and the
         * segments joining them.
         * @param points_name Polyscope name for the point cloud.
         * @param net_name Polyscope name for the segments.
         * @param points The control points.
         * @param segments Index pairs joining them.
         * @param color Color of both.
         * @param visible Whether the pair is currently shown.
         */
        void refresh_control_net(const char *points_name,
                                 const char *net_name,
                                 const std::vector<std::array<double, 3>> &points,
                                 const std::vector<std::array<int, 2>> &segments,
                                 const std::array<float, 3> &color,
                                 bool visible);
        /** @brief Creates a hex block spanning the model's bounding box, scaled by @p AMargin. */
        void create_bounding_box(double margin);
        /** @brief Registers (or removes) the mesh the blocking generates, at `m_subdivisions`. */
        void refresh_mesh();
        /** @brief Draws the "BIY operations" window: the parts of Polyscope's own panel worth
         * keeping, then biy's controls. */
        void draw_operations_panel();
        /** @brief Draws the "Scene" window, which replaces Polyscope's "Structures" panel with a
         * Geometry / Blocking / Mesh tree. */
        void draw_scene_panel();
        /**
         * @brief Draws the orientation gizmo: a small always-on-top compass, bottom-right of the 3D
         * view, showing the 6 world axes projected through the current camera frame — so the scene's
         * orientation reads at a glance, however the view has been rotated — with each axis dot
         * clickable to fly the camera to look along it.
         */
        void draw_gizmo();
        /** @brief Flies the camera to look along @p axis_dir at the scene's center, at Polyscope's
         * own home-view framing distance. */
        static void fly_to_axis_view(const glm::vec3 &axis_dir);
        /**
         * @brief Draws one Scene subsection: a named row delegating to the structure's own
         * Polyscope UI, so the usual visibility controls come from Polyscope itself.
         * @param label The subsection's name.
         * @param structure The structure to expose, or nullptr when the model has none of that kind
         *        (a boundary-representation model has no volume, say), which draws a greyed row.
         */
        static void draw_scene_entry(const char *label, polyscope::Structure *structure);
        /** @brief Draws the button panel. Assumes the ImGui frame is already set up. */
        void draw_panel();
        /** @brief Starts/continues/ends a corner drag from the current mouse state. */
        void handle_drag();
        /** @brief Tracks what the cursor is over in Cut mode, previewing the sheet it would cut, and
         * performs the cut on a click or on space. */
        void handle_cut();
        /** @brief Cuts the sheet the preview is currently showing, reporting what happened both in
         * the status line and on the terminal.
         * @param trigger How the cut was asked for, for the terminal line ("click" or "space"). */
        void perform_cut(const char *trigger);
        /** @brief Takes one edit back, or puts one back, and clears everything that described the
         * blocking as it was.
         * @param ARedo true to redo, false to undo. */
        void step_history(bool ARedo);

        /** @brief Tracks what the cursor is over in Collapse mode, lighting the sheet it would take
         * out, and collapses it on a click or on space. */
        void handle_collapse();
        /** @brief Collapses the sheet the preview is currently showing, reporting what happened both
         * in the status line and on the terminal.
         * @param trigger How the collapse was asked for, for the terminal line. */
        void perform_collapse(const char *trigger);
        /** @brief Tracks which block the cursor is over in Delete mode, highlighting it, and deletes
         * it on a click or on space. */
        void handle_delete();

        /** @brief Gathers block faces into a nappe in Pillow mode — a click adds one or takes it
         * back — and inserts the layer along it on space. */
        void handle_pillow();
        /** @brief Inserts a layer of blocks along the nappe gathered so far, reporting what happened
         * both in the status line and on the terminal.
         * @param trigger How it was asked for, for the terminal line. */
        void perform_pillow(const char *trigger);
        /** @brief Tracks which face the cursor is over in Chord mode, and which of its corners the
         * fold would run through, folding its chord away on a click or on space. */
        void handle_chord();
        /** @brief Folds the chord the preview is currently showing.
         * @param trigger How it was asked for, for the terminal line. */
        void perform_chord(const char *trigger);
        /** @brief Collects an edge and then the 2 faces its fan is cut at, in Open mode, and opens
         * the chord along them as soon as the second is picked. */
        void handle_open();
        /** @brief Opens the chord along what has been picked.
         * @param trigger How it was asked for, for the terminal line. */
        void perform_open(const char *trigger);
        /** @brief Re-reads which face the cursor is over, and which of its corners is nearest the
         * point picked — the 2 things a fold is named by.
         * @param screen_coords Current mouse position.
         * @return The face's id and that corner's node id, or nothing when the cursor is over no
         *         block face. */
        std::optional<std::pair<int, int>> pick_face(glm::vec2 screen_coords);
        /** @brief Which block edge a point in space is nearest, and how far along it.
         * @param target The point, as Polyscope reports it for a pick.
         * @return That edge's *position* in the display order and the parameter along it, or nothing
         *         when the blocking has no edge to measure against. */
        std::optional<std::pair<int, double>> nearest_edge(const glm::vec3 &target) const;
        /**
         * @brief Colours whichever faces are currently spoken for — the nappe being gathered, the
         * face a fold is aimed at, or the cuts an opening has so far.
         *
         * Painted onto the faces themselves rather than drawn as a second surface over them. A
         * highlight laid on top would sit in exactly the same plane as what it marks, and 2 coplanar
         * surfaces fight over the depth buffer: the mark shows in patches or not at all, which is
         * what the first version of this did. Colouring the face is also what "selected" is expected
         * to look like.
         */
        void refresh_face_preview();
        /** @brief The mean of one block's 8 corners.
         * @param block_id The block.
         * @return Its centre, or nothing when no block carries that id. */
        std::optional<glm::vec3> block_centre(int block_id);
        /** @brief Forgets every face, edge and corner picked for an operation in progress. */
        void clear_selection();
        /** @brief Re-reads which block the cursor is over.
         * @param screen_coords Current mouse position. */
        void update_delete_hover(glm::vec2 screen_coords);
        /** @brief Registers (or removes, when nothing is hovered) the highlight on the block about to
         * be deleted. */
        void refresh_delete_preview();
        /** @brief Deletes the block the highlight is currently showing, reporting what happened both
         * in the status line and on the terminal.
         * @param trigger How the deletion was asked for, for the terminal line. */
        void perform_delete(const char *trigger);
        /** @brief Re-reads which edge the cursor is over and rebuilds the cut preview from it.
         * @param screen_coords Current mouse position. */
        void update_cut_hover(glm::vec2 screen_coords);
        /** @brief Registers (or removes, when nothing is hovered) the sheet highlight and the markers
         * showing where the cut would land. */
        void refresh_cut_preview();
        /** @brief Where the cell of dimension @p ADim carrying @p AId sits in that dimension's
         * display arrays, or -1 if it is gone.
         * @param ADim 1 for edges, 3 for blocks.
         * @param AId The cell id.
         * @return Its position, for indexing the display arrays. */
        int display_position(int ADim, int AId) const;

        /** @brief A visually distinct colour for each block, from its index alone.
         * @param index The block's position in the blocking's own block order.
         * @return Its colour. */
        static glm::vec3 block_color(int index);
        /** @brief Reprojects screen coordinates onto the camera-facing plane through @p AAnchor. */
        static glm::vec3 screen_to_plane(glm::vec2 screen_coords, const glm::vec3 &anchor);

        std::mutex m_mutex;
        BiyConfig m_config;
        std::unique_ptr<app::GeomModelFacade> m_model;
        std::unique_ptr<app::BlockingFacade> m_blocking;
        std::string m_model_path;
        /** @brief Set by `request_refresh()`, consumed by the next `per_frame()`. Atomic rather than
         * guarded by `m_mutex`: the console calls `request_refresh()` while already holding that
         * mutex, and re-locking a non-recursive `std::mutex` would deadlock. */
        std::atomic<bool> m_view_dirty{false};
        /** @brief What the left mouse button currently does. */
        MouseMode m_mode = MouseMode::Camera;
        /** @brief Node id of the corner currently being dragged, if any. */
        std::optional<int> m_dragged_node;
        /** @brief Edge the cursor is currently over in Cut or Collapse mode, as its *id*, or unset
         * when the cursor is over no edge.
         *
         * An id rather than the position the pick came back as, because it is held across frames and
         * acted on later: a position stops meaning this edge the moment an operation renumbers the
         * traversal, and the operations here are what does that. What gets *drawn* is still
         * positional — see `m_sheet`. */
        std::optional<int> m_hover_edge;
        /** @brief The sheet that edge belongs to, as *positions* in the display arrays — this one is
         * drawn rather than acted on, and a renderer indexes a flat array. Empty when the sheet
         * cannot be cut homogeneously, which is itself worth showing rather than hiding. */
        std::vector<int> m_sheet;
        /** @brief Where along the hovered edge the cut would fall. Driven by the cursor, and editable
         * as a number in the panel for a cut that has to land on an exact value. */
        float m_cut_param = 0.5f;
        /** @brief Last position the cut hover was computed for, so a still cursor costs nothing:
         * each hover test is a full Polyscope pick, which is a render pass. */
        glm::vec2 m_last_cut_mouse{-1.0f, -1.0f};
        /** @brief Block the cursor is currently over in Delete mode, as its *id*, or unset when the
         * cursor is over none. @see m_hover_edge for why an id and not a position. */
        std::optional<int> m_hover_block;
        /** @copydoc m_last_cut_mouse */
        glm::vec2 m_last_delete_mouse{-1.0f, -1.0f};
        /** @brief The nappe gathered so far in Pillow mode, as face *ids*: it is held across frames
         * and acted on later. @see m_hover_edge */
        std::vector<int> m_nappe;
        /** @brief Face the cursor is over in Chord mode, as its id, or unset when over none. */
        std::optional<int> m_hover_face;
        /** @brief Which corner of that face the fold would run through — the one nearest the point
         * picked, so aiming at a corner of a face is what names the diagonal. */
        std::optional<int> m_hinge_node;
        /** @brief Edge picked in Open mode, as its id, before its 2 cut faces are. */
        std::optional<int> m_open_edge;
        /** @brief The faces picked so far to cut that edge's fan at, as ids; 2 of them and the
         * opening runs. */
        std::vector<int> m_open_faces;
        /** @copydoc m_last_cut_mouse */
        glm::vec2 m_last_face_mouse{-1.0f, -1.0f};
        /** @brief Whether what is marked on the block faces still matches what is selected.
         *
         * The mark is repainted once a frame rather than at each place that changes it, so that a
         * rebuild of the Polyscope structures — which drops the mark with them — needs no one to
         * remember to put it back. This is what stops that costing a full pass over the face grid on
         * every frame in which nothing happened. */
        bool m_selection_dirty = true;
        /** @brief How far a pillow's inside, or an opening's 2 halves, move — as a fraction of the
         * mean length of the edges at each corner that moves, so a layer follows the size of the
         * blocks it is inserted between. */
        float m_thickness = 0.25f;
        /** @brief Subdivisions used when displaying the blocking (1 = raw block structure). */
        int m_subdivisions = 1;
        /** @brief Whether the block structure's own edges are currently drawn. */
        bool m_show_block_edges = true;
        /** @brief How the generated mesh is currently coloured. */
        MeshColoring m_mesh_coloring = MeshColoring::Uniform;
        /** @brief Whether the generated mesh is currently drawn. */
        bool m_show_mesh = false;
        /** @brief Whether each kind of control display is currently drawn. */
        bool m_show_edge_control = false;
        /** @copydoc m_show_edge_control */
        bool m_show_face_control = false;
        /** @copydoc m_show_edge_control */
        bool m_show_block_control = false;
        /** @brief Samples per block edge; > 1 traces a curved edge rather than its chord. */
        int m_edge_samples = 8;
        /** @brief Edge order the blocking is currently built at. Changeable from the panel — the
         * degree is carried by the geometry rather than by its C++ type, so raising or lowering it
         * refits the structure in place instead of needing a new one. */
        int m_order = 3;
        /** @brief Per-dimension snapping tolerances, shared by the "Classify" button and by the
         * snap that runs when a dragged corner is released. Separate values because the scales
         * differ: one loose enough to catch a surface would snap corners to the wrong vertex. */
        float m_tol_vertex = 0.1f;
        /** @copydoc m_tol_vertex */
        float m_tol_curve = 0.1f;
        /** @copydoc m_tol_vertex */
        float m_tol_surface = 0.1f;
        /** @brief Last status/error line shown in the panel. */
        std::string m_status;
    };

} // namespace gecko::biy
