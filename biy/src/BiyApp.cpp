#include "BiyApp.h"

#include <algorithm>
#include <iostream>
#include <limits>

#include <imgui.h>
#include <polyscope/options.h>
#include <polyscope/pick.h>
#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/view.h>
#include <polyscope/volume_mesh.h>

namespace gecko::biy {

    namespace {
        /** @brief Polyscope structure names, also used to recognize what a pick hit. */
        constexpr const char *MODEL_SURFACE = "model surface";
        constexpr const char *MODEL_VOLUME = "model volume";
        constexpr const char *BLOCK_CORNERS = "block corners";
        constexpr const char *BLOCK_QUADS = "block quads";
        constexpr const char *BLOCK_HEXES = "block hexes";
        /** @brief Overlay holding just the corner being dragged, drawn bigger and in its own color.
         * A separate structure rather than a per-point quantity on BLOCK_CORNERS: it keeps the
         * highlight's radius and color independent of the base cloud's, and picking still reports
         * the base cloud underneath. */
        constexpr const char *DRAG_HIGHLIGHT = "dragged corner";

        glm::vec3 to_glm(const std::array<float, 3> &c) { return {c[0], c[1], c[2]}; }
    } // namespace

    BiyApp::BiyApp(const std::string &model_path)
        : m_model(std::make_unique<python::GeomModelFacade>(model_path)), m_model_path(model_path) {
        std::string config_message;
        m_config = BiyConfig::load("biy_config.json", config_message);
        std::cout << config_message << "\n";
        register_model();
        m_status = "Loaded " + model_path + " — " + config_message;
    }

    void BiyApp::set_mouse_mode(MouseMode mode) {
        m_mode = mode;
        // Switched here, a whole frame before any drag is processed — Polyscope reads this at the
        // top of its frame, before the user callback gets a say.
        polyscope::options::doDefaultMouseInteraction = (mode == MouseMode::Camera);
        if (mode == MouseMode::Camera && m_dragged_node) {
            m_dragged_node.reset();
            show_highlight(std::nullopt);
        }
    }

    void BiyApp::register_model() {
        const auto vertices = m_model->mesh_vertices();
        const auto triangles = m_model->mesh_triangles();
        const auto tets = m_model->mesh_tets();

        if (!triangles.empty()) {
            polyscope::registerSurfaceMesh(MODEL_SURFACE, vertices, triangles);
        }
        if (!tets.empty()) {
            polyscope::registerTetMesh(MODEL_VOLUME, vertices, tets);
            // A tet model's own boundary triangles already show its shape; keeping the volume
            // visible on top of them just hides the blocking being built inside it.
            polyscope::getVolumeMesh(MODEL_VOLUME)->setEnabled(false);
        }
        polyscope::view::resetCameraToHomeView();
    }

    python::BlockingFacade &BiyApp::blocking() {
        if (!m_blocking) m_blocking = std::make_unique<python::BlockingFacade>(*m_model);
        return *m_blocking;
    }

    void BiyApp::refresh_view() {
        if (!m_blocking) return;

        const auto node_ids = m_blocking->node_ids();
        std::vector<glm::vec3> corners;
        corners.reserve(node_ids.size());
        for (const int id : node_ids) {
            const auto p = m_blocking->node_position(id);
            corners.emplace_back(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]));
        }
        auto *cloud = polyscope::registerPointCloud(BLOCK_CORNERS, corners);
        cloud->setPointRadius(m_config.corner_radius);
        cloud->setPointColor(to_glm(m_config.corner_color));
        if (m_dragged_node) show_highlight(m_dragged_node);

        const auto vertices = m_blocking->mesh_vertices(m_subdivisions);
        const auto quads = m_blocking->mesh_quads(m_subdivisions);
        const auto hexes = m_blocking->mesh_hexes(m_subdivisions);

        if (!quads.empty()) {
            polyscope::registerSurfaceMesh(BLOCK_QUADS, vertices, quads);
        } else if (polyscope::hasSurfaceMesh(BLOCK_QUADS)) {
            polyscope::removeStructure(polyscope::getSurfaceMesh(BLOCK_QUADS));
        }

        if (!hexes.empty()) {
            polyscope::registerHexMesh(BLOCK_HEXES, vertices, hexes);
        } else if (polyscope::hasVolumeMesh(BLOCK_HEXES)) {
            polyscope::removeStructure(polyscope::getVolumeMesh(BLOCK_HEXES));
        }
    }

    void BiyApp::show_highlight(std::optional<int> node_id) {
        if (!node_id || !m_blocking) {
            if (polyscope::hasPointCloud(DRAG_HIGHLIGHT)) {
                polyscope::removeStructure(polyscope::getPointCloud(DRAG_HIGHLIGHT));
            }
            return;
        }

        const auto p = m_blocking->node_position(*node_id);
        const std::vector<glm::vec3> single{
            {static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2])}};
        auto *cloud = polyscope::registerPointCloud(DRAG_HIGHLIGHT, single);
        cloud->setPointRadius(m_config.corner_highlight_radius);
        cloud->setPointColor(to_glm(m_config.corner_highlight_color));
    }

    void BiyApp::create_bounding_box(double margin) {
        const auto vertices = m_model->mesh_vertices();
        if (vertices.empty()) {
            m_status = "Model has no vertices to bound";
            return;
        }

        std::array<double, 3> lo = vertices.front();
        std::array<double, 3> hi = vertices.front();
        for (const auto &v : vertices) {
            for (std::size_t k = 0; k < 3; ++k) {
                lo[k] = std::min(lo[k], v[k]);
                hi[k] = std::max(hi[k], v[k]);
            }
        }
        for (std::size_t k = 0; k < 3; ++k) {
            const double mid = 0.5 * (lo[k] + hi[k]);
            const double half = 0.5 * (hi[k] - lo[k]) * (1.0 + margin);
            lo[k] = mid - half;
            hi[k] = mid + half;
        }

        // HEX8 ordering: bottom perimeter 0-1-2-3, top perimeter 4-5-6-7 directly above.
        const std::vector<std::array<double, 3>> corners{{lo[0], lo[1], lo[2]},
                                                          {hi[0], lo[1], lo[2]},
                                                          {hi[0], hi[1], lo[2]},
                                                          {lo[0], hi[1], lo[2]},
                                                          {lo[0], lo[1], hi[2]},
                                                          {hi[0], lo[1], hi[2]},
                                                          {hi[0], hi[1], hi[2]},
                                                          {lo[0], hi[1], hi[2]}};
        blocking().create_hex_block(corners);
        refresh_view();
        m_status = "Created bounding box";
    }

    void BiyApp::per_frame() {
        const std::lock_guard<std::mutex> lock(m_mutex);
        draw_panel();
        handle_drag();
    }

    void BiyApp::draw_panel() {
        ImGui::TextUnformatted(m_model_path.c_str());
        ImGui::Separator();

        // Mode: keyboard shortcuts first, so they work wherever the mouse happens to be — but not
        // while a widget has the keyboard, or typing "1e-3" into a tolerance field would switch mode
        // on the "e".
        if (!ImGui::GetIO().WantCaptureKeyboard) {
            if (ImGui::IsKeyPressed(ImGuiKey_C)) set_mouse_mode(MouseMode::Camera);
            if (ImGui::IsKeyPressed(ImGuiKey_E)) set_mouse_mode(MouseMode::Edit);
        }

        ImGui::TextUnformatted("Mouse mode");
        if (ImGui::RadioButton("Camera (C)", m_mode == MouseMode::Camera)) set_mouse_mode(MouseMode::Camera);
        ImGui::SameLine();
        if (ImGui::RadioButton("Edit (E)", m_mode == MouseMode::Edit)) set_mouse_mode(MouseMode::Edit);
        ImGui::TextWrapped(m_mode == MouseMode::Edit
                               ? "Drag a block corner to move it. Camera navigation is off."
                               : "Rotate/pan/zoom the view. Switch to Edit to move corners.");
        ImGui::Separator();

        if (ImGui::Button("Create bounding box")) create_bounding_box(0.1);

        const bool has_blocking = m_blocking != nullptr;
        ImGui::BeginDisabled(!has_blocking);

        if (ImGui::Button("Build connectivity")) {
            m_blocking->build_connectivity();
            refresh_view();
            m_status = "Built connectivity";
        }

        ImGui::SetNextItemWidth(120.0f);
        ImGui::InputFloat("classify tol", &m_classify_tol);
        if (ImGui::Button("Classify")) {
            m_blocking->classify(m_classify_tol);
            refresh_view();
            m_status = "Classified onto the model";
        }

        ImGui::SetNextItemWidth(120.0f);
        if (ImGui::InputInt("subdivisions", &m_subdivisions)) {
            m_subdivisions = std::max(1, m_subdivisions);
            refresh_view();
        }

        if (ImGui::Button("Export VTK")) {
            const std::string out = "biy_blocking.vtk";
            m_blocking->write_vtk(m_subdivisions, out);
            m_status = "Wrote " + out;
        }

        if (has_blocking) {
            ImGui::Text("nodes %zu  edges %zu  faces %zu  blocks %zu",
                        m_blocking->nb_cells(0),
                        m_blocking->nb_cells(1),
                        m_blocking->nb_cells(2),
                        m_blocking->nb_cells(3));
        }
        ImGui::EndDisabled();

        ImGui::Separator();
        ImGui::TextWrapped("Drag a block corner with the left mouse button. "
                           "The Python console on stdin drives the same `model` and `blocking`.");
        ImGui::TextUnformatted(m_status.c_str());
    }

    glm::vec3 BiyApp::screen_to_plane(glm::vec2 screen_coords, const glm::vec3 &anchor) {
        const glm::vec3 eye = polyscope::view::getCameraWorldPosition();
        const glm::vec3 ray = polyscope::view::screenCoordsToWorldRay(screen_coords);
        // Plane through `anchor`, facing the camera: its normal is the view direction, so the drag
        // tracks the mouse exactly while keeping the corner at its original depth.
        const glm::vec3 normal = glm::normalize(eye - anchor);
        const float denom = glm::dot(ray, normal);
        if (std::abs(denom) < 1e-6f) return anchor; // ray parallel to the plane: nothing sensible to do
        const float t = glm::dot(anchor - eye, normal) / denom;
        return eye + ray * t;
    }

    void BiyApp::handle_drag() {
        if (!m_blocking || m_mode != MouseMode::Edit) return;
        ImGuiIO &io = ImGui::GetIO();

        if (!m_dragged_node && io.MouseClicked[0] && !io.WantCaptureMouse) {
            const glm::vec2 screen_coords{io.MousePos.x, io.MousePos.y};
            const polyscope::PickResult pick = polyscope::pickAtScreenCoords(screen_coords);
            if (pick.isHit && pick.structureName == BLOCK_CORNERS) {
                const auto node_ids = m_blocking->node_ids();
                if (pick.localIndex < node_ids.size()) {
                    m_dragged_node = node_ids[pick.localIndex];
                    show_highlight(m_dragged_node);
                    m_status = "Dragging corner " + std::to_string(*m_dragged_node);
                }
            }
        }

        if (m_dragged_node) {
            if (io.MouseDown[0]) {
                const auto current = m_blocking->node_position(*m_dragged_node);
                const glm::vec3 anchor(static_cast<float>(current[0]),
                                       static_cast<float>(current[1]),
                                       static_cast<float>(current[2]));
                const glm::vec3 target = screen_to_plane({io.MousePos.x, io.MousePos.y}, anchor);
                m_blocking->move_node(*m_dragged_node, target.x, target.y, target.z);
                refresh_view();
            } else {
                m_status = "Moved corner " + std::to_string(*m_dragged_node);
                m_dragged_node.reset();
                show_highlight(std::nullopt);
            }
        }
    }

} // namespace gecko::biy
