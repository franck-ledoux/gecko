#pragma once

#include <memory>
#include <mutex>
#include <optional>
#include <string>

#include <glm/glm.hpp>

#include "BiyConfig.h"
#include "BlockingFacade.h"
#include "GeomModelFacade.h"

namespace gecko::biy {

    /**
     * @brief What the left mouse button does. Polyscope processes its own camera navigation
     * *before* the per-frame user callback runs, so a drag can't be intercepted after the fact —
     * navigation has to be switched off ahead of time, which makes this a genuine mode rather than
     * a per-event decision.
     */
    enum class MouseMode {
        Camera, ///< Polyscope's usual navigation: rotate/pan/zoom the view.
        Edit    ///< Navigation off; dragging picks up and moves a block corner.
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
        /**
         * @brief Constructor. Loads the geometric model and registers it for display.
         * @param model_path Path to the .msh file to load.
         * @throw std::runtime_error if the file cannot be read.
         */
        explicit BiyApp(const std::string &model_path);

        /** @brief Draws this frame's ImGui panel and services any in-progress corner drag. Called
         * by Polyscope once per frame, on the main thread. */
        void per_frame();

        /** @brief The lock guarding the facades, shared with the Python console. @return The mutex. */
        [[nodiscard]] std::mutex &mutex() { return m_mutex; }

        /** @brief The geometric model, for injection into the Python console. @return The model. */
        [[nodiscard]] python::GeomModelFacade &model() { return *m_model; }

        /**
         * @brief Gets the blocking, creating an empty one on first use so the Python console always
         * has a `blocking` to talk to.
         * @return The blocking.
         */
        [[nodiscard]] python::BlockingFacade &blocking();

        /** @brief Re-reads the facades and refreshes every Polyscope structure. Call after any
         * mutation — including ones made from the Python console. */
        void refresh_view();

    private:
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
        /** @brief Creates a hex block spanning the model's bounding box, scaled by @p AMargin. */
        void create_bounding_box(double margin);
        /** @brief Draws the button panel. Assumes the ImGui frame is already set up. */
        void draw_panel();
        /** @brief Starts/continues/ends a corner drag from the current mouse state. */
        void handle_drag();
        /** @brief Reprojects screen coordinates onto the camera-facing plane through @p AAnchor. */
        static glm::vec3 screen_to_plane(glm::vec2 screen_coords, const glm::vec3 &anchor);

        std::mutex m_mutex;
        BiyConfig m_config;
        std::unique_ptr<python::GeomModelFacade> m_model;
        std::unique_ptr<python::BlockingFacade> m_blocking;
        std::string m_model_path;
        /** @brief What the left mouse button currently does. */
        MouseMode m_mode = MouseMode::Camera;
        /** @brief Node id of the corner currently being dragged, if any. */
        std::optional<int> m_dragged_node;
        /** @brief Subdivisions used when displaying the blocking (1 = raw block structure). */
        int m_subdivisions = 1;
        /** @brief Tolerance the "Classify" button passes to Blocking::classify(). */
        float m_classify_tol = 0.1f;
        /** @brief Last status/error line shown in the panel. */
        std::string m_status;
    };

} // namespace gecko::biy
