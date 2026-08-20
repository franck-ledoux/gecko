#include "BiyApp.h"

#include <algorithm>
#include <iostream>
#include <limits>

#include <imgui.h>
#include <polyscope/curve_network.h>
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
        /** @brief The block faces themselves, sampled into quads — not the generated mesh's quads,
         * which exist only for standalone 2D blocks (see `refresh_view()`). */
        constexpr const char *BLOCK_FACES = "block faces";
        constexpr const char *BLOCK_HEXES = "block hexes";
        /** @brief The control points driving curved edges, and the polygon joining them. */
        constexpr const char *CONTROL_POINTS = "control points";
        constexpr const char *CONTROL_POLYGON = "control polygon";
        /** @brief Overlay holding just the corner being dragged, drawn bigger and in its own color.
         * A separate structure rather than a per-point quantity on BLOCK_CORNERS: it keeps the
         * highlight's radius and color independent of the base cloud's, and picking still reports
         * the base cloud underneath. */
        constexpr const char *DRAG_HIGHLIGHT = "dragged corner";
        /** @brief The block structure's own edges, traced along their curves — distinct from the
         * subdivision lines of the mesh BLOCK_QUADS/BLOCK_HEXES show. */
        constexpr const char *BLOCK_EDGES = "block edges";

        glm::vec3 to_glm(const std::array<float, 3> &c) { return {c[0], c[1], c[2]}; }

        /** @brief What a classification dimension means, for the status line. */
        const char *classification_name(int classification_dim) {
            switch (classification_dim) {
                case 0:
                    return "a vertex";
                case 1:
                    return "a curve";
                case 2:
                    return "a surface";
                case 3:
                    return "a volume";
                default:
                    return "nothing";
            }
        }
    } // namespace

    BiyApp::BiyApp(const std::string &model_path, int order)
        : m_model(std::make_unique<python::GeomModelFacade>(model_path)), m_model_path(model_path), m_order(order) {
        std::string config_message;
        m_config = BiyConfig::load("biy_config.json", config_message);
        m_show_block_edges = m_config.show_block_edges;
        m_show_control_points = m_config.show_control_points;
        m_tol_vertex = static_cast<float>(m_config.tol_vertex);
        m_tol_curve = static_cast<float>(m_config.tol_curve);
        m_tol_surface = static_cast<float>(m_config.tol_surface);
        std::cout << config_message << "\n";
        register_model();
        m_status = "Loaded " + model_path + " (order " + std::to_string(order) + ") — " + config_message;
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
            auto *surface = polyscope::registerSurfaceMesh(MODEL_SURFACE, vertices, triangles);
            // See-through by default: a blocking is built around and inside its model, and an
            // opaque model hides the very thing being edited.
            surface->setTransparency(static_cast<float>(m_config.model_transparency));
        }
        if (!tets.empty()) {
            polyscope::registerTetMesh(MODEL_VOLUME, vertices, tets);
            // A tet model's own boundary triangles already show its shape; keeping the volume
            // visible on top of them just hides the blocking being built inside it.
            polyscope::getVolumeMesh(MODEL_VOLUME)->setEnabled(false);
        }
        polyscope::view::resetCameraToHomeView();
        freeze_scene_extents();
    }

    void BiyApp::freeze_scene_extents() {
        // Called once the model is registered, so the scene's bounding box and length scale are the
        // model's. From here on Polyscope stops recomputing them as structures come and go — which
        // is what keeps the ground plane still: its height is derived from that bounding box, so
        // without this, creating a block or dragging one corner outwards slides the ground (and
        // rescales the view) out from under the model. The blocking is meant to be edited freely;
        // the scene it's edited in should not move in response.
        polyscope::options::automaticallyComputeSceneExtents = false;
    }

    python::BlockingFacade &BiyApp::blocking() {
        if (!m_blocking) m_blocking = std::make_unique<python::BlockingFacade>(*m_model, m_order);
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

        // One color per corner, by what it is classified on — the state of a blocking being fitted
        // to a model is mostly "which corners have found their geometry yet", so it belongs on the
        // corners themselves rather than in a side panel.
        const auto dims = m_blocking->node_classification_dims();
        std::vector<glm::vec3> colors;
        colors.reserve(dims.size());
        for (const int dim : dims)
            colors.push_back(to_glm(m_config.color_for(dim)));
        cloud->addColorQuantity("classification", colors)->setEnabled(true);

        if (m_dragged_node) show_highlight(m_dragged_node);

        const auto edge_points = m_blocking->edge_vertices(m_edge_samples);
        if (!edge_points.empty()) {
            auto *edges =
                polyscope::registerCurveNetwork(BLOCK_EDGES, edge_points, m_blocking->edge_segments(m_edge_samples));
            edges->setRadius(m_config.block_edge_radius);
            edges->setColor(to_glm(m_config.block_edge_color));
            edges->setEnabled(m_show_block_edges);
            // One color per edge, spread over the `m_edge_samples` segments it was sampled into.
            apply_classification_colors(edges, m_blocking->edge_classification_dims(), [this](auto *s, const auto &c) {
                std::vector<glm::vec3> per_segment;
                per_segment.reserve(c.size() * static_cast<std::size_t>(m_edge_samples));
                for (const auto &color : c) {
                    per_segment.insert(per_segment.end(), static_cast<std::size_t>(m_edge_samples), color);
                }
                return s->addEdgeColorQuantity("classification", per_segment);
            });
        } else if (polyscope::hasCurveNetwork(BLOCK_EDGES)) {
            polyscope::removeStructure(polyscope::getCurveNetwork(BLOCK_EDGES));
        }

        refresh_control_points();

        // Block faces get their own structure rather than riding on the generated mesh: to_mesh()
        // only emits quads for standalone 2D blocks, so a hex block's 6 bounding faces would
        // otherwise be invisible — and unclassifiable to the eye.
        const auto face_points = m_blocking->face_grid_vertices(m_subdivisions);
        if (!face_points.empty()) {
            auto *faces =
                polyscope::registerSurfaceMesh(BLOCK_FACES, face_points, m_blocking->face_grid_quads(m_subdivisions));
            const auto owners = m_blocking->face_grid_owners(m_subdivisions);
            apply_classification_colors(
                faces, m_blocking->face_classification_dims(), [&owners](auto *s, const auto &c) {
                    std::vector<glm::vec3> per_quad;
                    per_quad.reserve(owners.size());
                    for (const int owner : owners) {
                        per_quad.push_back(c[static_cast<std::size_t>(owner)]);
                    }
                    return s->addFaceColorQuantity("classification", per_quad);
                });
        } else if (polyscope::hasSurfaceMesh(BLOCK_FACES)) {
            polyscope::removeStructure(polyscope::getSurfaceMesh(BLOCK_FACES));
        }

        const auto hexes = m_blocking->mesh_hexes(m_subdivisions);
        if (!hexes.empty()) {
            polyscope::registerHexMesh(BLOCK_HEXES, m_blocking->mesh_vertices(m_subdivisions), hexes);
        } else if (polyscope::hasVolumeMesh(BLOCK_HEXES)) {
            polyscope::removeStructure(polyscope::getVolumeMesh(BLOCK_HEXES));
        }
    }

    template<typename TStructure, typename TAdd>
    void BiyApp::apply_classification_colors(TStructure *structure, const std::vector<int> &dims, TAdd add) {
        std::vector<glm::vec3> colors;
        colors.reserve(dims.size());
        for (const int dim : dims) {
            colors.push_back(to_glm(m_config.color_for(dim)));
        }
        add(structure, colors)->setEnabled(true);
    }

    void BiyApp::refresh_control_points() {
        // A degree-1 edge's 2 control points are just its endpoints, already drawn as corners —
        // there is nothing to reveal, so the structure is never created at that order.
        const auto points =
            (m_blocking && m_order > 1) ? m_blocking->edge_control_points() : std::vector<std::array<double, 3>>{};
        if (points.empty()) {
            if (polyscope::hasPointCloud(CONTROL_POINTS)) {
                polyscope::removeStructure(polyscope::getPointCloud(CONTROL_POINTS));
            }
            if (polyscope::hasCurveNetwork(CONTROL_POLYGON)) {
                polyscope::removeStructure(polyscope::getCurveNetwork(CONTROL_POLYGON));
            }
            return;
        }

        auto *cloud = polyscope::registerPointCloud(CONTROL_POINTS, points);
        cloud->setPointRadius(m_config.control_point_radius);
        cloud->setPointColor(to_glm(m_config.control_point_color));
        cloud->setEnabled(m_show_control_points);

        auto *polygon = polyscope::registerCurveNetwork(CONTROL_POLYGON, points, m_blocking->edge_control_polygons());
        polygon->setRadius(m_config.control_polygon_radius);
        polygon->setColor(to_glm(m_config.control_point_color));
        polygon->setEnabled(m_show_control_points);
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

    void BiyApp::request_refresh() { m_view_dirty = true; }

    void BiyApp::per_frame() {
        const std::lock_guard<std::mutex> lock(m_mutex);
        // Any refresh asked for off this thread happens here, before anything is drawn: Polyscope
        // renders outside our callback and outside this lock, so registering structures from the
        // console thread would race with the very loop drawing them.
        if (m_view_dirty.exchange(false)) {
            refresh_view();
        }
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
        ImGui::TextWrapped(m_mode == MouseMode::Edit ? "Drag a block corner to move it. Camera navigation is off."
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

        // Snapping tolerances, shared by the Classify button and by the snap on drag release.
        ImGui::SetNextItemWidth(90.0f);
        ImGui::InputFloat("tol vertex", &m_tol_vertex);
        ImGui::SetNextItemWidth(90.0f);
        ImGui::InputFloat("tol curve", &m_tol_curve);
        ImGui::SetNextItemWidth(90.0f);
        ImGui::InputFloat("tol surface", &m_tol_surface);
        if (ImGui::Button("Classify")) {
            m_blocking->classify(m_tol_vertex, m_tol_curve, m_tol_surface);
            refresh_view();
            m_status = "Classified onto the model";
        }

        ImGui::SetNextItemWidth(120.0f);
        if (ImGui::InputInt("subdivisions", &m_subdivisions)) {
            m_subdivisions = std::max(1, m_subdivisions);
            refresh_view();
        }

        if (ImGui::Checkbox("Show block edges", &m_show_block_edges)) refresh_view();

        ImGui::BeginDisabled(m_order <= 1);
        if (ImGui::Checkbox("Show control points", &m_show_control_points)) refresh_view();
        ImGui::EndDisabled();
        if (m_order <= 1 && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
            ImGui::SetTooltip("Order 1 edges are straight: their control points are just the corners.");
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

            ImGui::TextUnformatted("Corners:");
            const auto swatch = [](const std::array<float, 3> &c, const char *label) {
                ImGui::SameLine();
                ImGui::TextColored(ImVec4(c[0], c[1], c[2], 1.0f), "%s", label);
            };
            swatch(m_config.corner_color_unclassified, "free");
            swatch(m_config.corner_color_on_vertex, "vertex");
            swatch(m_config.corner_color_on_curve, "curve");
            swatch(m_config.corner_color_on_surface, "surface");
            swatch(m_config.corner_color_on_volume, "volume");
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
                const glm::vec3 anchor(
                    static_cast<float>(current[0]), static_cast<float>(current[1]), static_cast<float>(current[2]));
                const glm::vec3 target = screen_to_plane({io.MousePos.x, io.MousePos.y}, anchor);
                m_blocking->move_node(*m_dragged_node, target.x, target.y, target.z);
                refresh_view();
            } else {
                // Snap on release: the corner settles onto whatever it landed near, and every
                // edge/face touching it is reclassified and refitted to match — so the colors on
                // screen tell the truth again the moment the button comes up.
                const int released = *m_dragged_node;
                m_blocking->snap_node(released, m_tol_vertex, m_tol_curve, m_tol_surface);
                m_dragged_node.reset();
                show_highlight(std::nullopt);
                refresh_view();

                const auto dims = m_blocking->node_classification_dims();
                const auto ids = m_blocking->node_ids();
                const auto it = std::find(ids.begin(), ids.end(), released);
                const int dim = (it != ids.end()) ? dims[static_cast<std::size_t>(it - ids.begin())] : -1;
                m_status = "Corner " + std::to_string(released) + " snapped onto " + classification_name(dim);
            }
        }
    }

} // namespace gecko::biy
