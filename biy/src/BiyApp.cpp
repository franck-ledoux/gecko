#include "BiyApp.h"

#include <algorithm>
#include <iostream>
#include <limits>

#include <imgui.h>
#include <polyscope/curve_network.h>
#include <polyscope/options.h>
#include <polyscope/internal.h>
#include <polyscope/pick.h>
#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>
#include <polyscope/render/engine.h>
#include <polyscope/screenshot.h>
#include <polyscope/structure.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/view.h>
#include <polyscope/volume_mesh.h>

namespace gecko::biy {

    namespace {
        // Polyscope structure names. They double as the labels shown in the Scene panel and as
        // what a pick reports, so they read as the section they belong to.

        /** @brief The geometric model, one structure per kind of entity — the Geometry section's 4
         * subsections. */
        constexpr const char *GEOM_VOLUMES = "volumes";
        constexpr const char *GEOM_SURFACES = "surfaces";
        constexpr const char *GEOM_CURVES = "curves";
        constexpr const char *GEOM_POINTS = "points";
        /** @brief The mesh the blocking generates, as opposed to the block structure itself. */
        constexpr const char *MESH_QUADS = "mesh quads";
        constexpr const char *MESH_HEXES = "mesh hexes";
        /** @brief The block structure itself — the Blocking section's 4 subsections. */
        constexpr const char *BLOCK_VERTICES = "vertices";
        /** @brief The block faces themselves, sampled into quads — not the generated mesh's quads,
         * which exist only for standalone 2D blocks (see `refresh_view()`). */
        constexpr const char *BLOCK_FACES = "faces";
        constexpr const char *BLOCK_BLOCKS = "blocks";
        /** @brief The control points driving each curved cell, and the polygon/net/lattice joining
         * them. Three separate pairs: a face's net already contains its edges' points and a block's
         * lattice contains both, so they are told apart by color and toggled independently. */
        constexpr const char *EDGE_CONTROL_POINTS = "edge control points";
        constexpr const char *EDGE_CONTROL_NET = "edge control polygon";
        constexpr const char *FACE_CONTROL_POINTS = "face control points";
        constexpr const char *FACE_CONTROL_NET = "face control net";
        constexpr const char *BLOCK_CONTROL_POINTS = "block control points";
        constexpr const char *BLOCK_CONTROL_NET = "block control lattice";
        /** @brief Overlay holding just the corner being dragged, drawn bigger and in its own color.
         * A separate structure rather than a per-point quantity on BLOCK_VERTICES: it keeps the
         * highlight's radius and color independent of the base cloud's, and picking still reports
         * the base cloud underneath. */
        constexpr const char *DRAG_HIGHLIGHT = "dragged corner";
        /** @brief The block structure's own edges, traced along their curves — distinct from the
         * subdivision lines of the generated mesh. */
        constexpr const char *BLOCK_EDGES = "edges";

        glm::vec3 to_glm(const std::array<float, 3> &c) { return {c[0], c[1], c[2]}; }

        // Polyscope looks structures up per type, so the Scene panel needs one accessor per kind.
        // Each returns null when that structure does not exist — a boundary-representation model has
        // no volume, a 2D blocking no blocks — which the panel shows as a greyed-out row.
        polyscope::Structure *volume_or_null(const char *name) {
            return polyscope::hasVolumeMesh(name) ? polyscope::getVolumeMesh(name) : nullptr;
        }
        polyscope::Structure *surface_or_null(const char *name) {
            return polyscope::hasSurfaceMesh(name) ? polyscope::getSurfaceMesh(name) : nullptr;
        }
        polyscope::Structure *curves_or_null(const char *name) {
            return polyscope::hasCurveNetwork(name) ? polyscope::getCurveNetwork(name) : nullptr;
        }
        polyscope::Structure *points_or_null(const char *name) {
            return polyscope::hasPointCloud(name) ? polyscope::getPointCloud(name) : nullptr;
        }

        /** @brief Whether @p structure is currently shown, or nullopt if it does not exist yet. */
        std::optional<bool> enabled_state(polyscope::Structure *structure) {
            return structure ? std::optional<bool>(structure->isEnabled()) : std::nullopt;
        }

        /**
         * @brief Restores a structure's visibility across a re-registration.
         *
         * refresh_view() rebuilds these structures on every edit, and registering over a name resets
         * it to visible — which would silently undo whatever was just ticked in the Scene panel. So
         * a previous state, when there is one, wins over the default.
         * @param structure The freshly registered structure.
         * @param previous What the same structure showed before, or nullopt if it is new.
         * @param fallback The state to use when it is new.
         */
        void restore_enabled(polyscope::Structure *structure, std::optional<bool> previous, bool fallback) {
            structure->setEnabled(previous.value_or(fallback));
        }

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
        m_show_mesh = m_config.show_mesh;
        m_show_edge_control = m_config.show_edge_control;
        m_show_face_control = m_config.show_face_control;
        m_show_block_control = m_config.show_block_control;
        m_tol_vertex = static_cast<float>(m_config.tol_vertex);
        m_tol_curve = static_cast<float>(m_config.tol_curve);
        m_tol_surface = static_cast<float>(m_config.tol_surface);
        std::cout << config_message << "\n";

        // biy draws its own panels in place of Polyscope's: no default windows, and no anonymous
        // wrapper around the per-frame callback, so the windows below can carry real titles.
        // leftWindowsWidth is normally initialised lazily by the very builders being dropped, so it
        // has to be set here or it stays at its -1 sentinel.
        polyscope::options::buildDefaultGuiPanels = false;
        polyscope::options::openImGuiWindowForUserCallback = false;
        polyscope::internal::leftWindowsWidth = static_cast<float>(m_config.panel_width) * polyscope::options::uiScale;
        // The right side is trickier: buildPickGui() (kept, for the "Selection" box) calls Polyscope's
        // own ensureWindowWidthsSet() every frame, which unconditionally *recomputes*
        // internal::rightWindowsWidth from options::rightGuiPaneWidth (500 by default) — so setting
        // internal::rightWindowsWidth directly, the way leftWindowsWidth is set above, would only
        // last until the next buildPickGui() call, and both Scene and Selection would drift to 500
        // regardless of panel_width. Setting the *source* option instead makes every later
        // recomputation land back on the width actually asked for.
        polyscope::options::rightGuiPaneWidth = static_cast<int>(m_config.panel_width);
        polyscope::internal::rightWindowsWidth = static_cast<float>(m_config.panel_width) * polyscope::options::uiScale;

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
        const auto curve_points = m_model->curve_vertices();
        const auto curve_segments = m_model->curve_segments();
        const auto points = m_model->vertex_positions();
        const auto volume_boundary = m_model->volume_boundary_triangles();

        // The model's 4 kinds of entity get 4 structures, matching the Geometry section's
        // subsections one for one — where it used to be a single undifferentiated triangle soup
        // plus a hidden tet mesh, with its curves and vertices not displayed at all.
        if (!volume_boundary.empty()) {
            // A genuine SurfaceMesh over the tet mesh's own boundary faces, not a VolumeMesh: unlike
            // what the name suggests, Polyscope's VolumeMesh always renders every tet face, interior
            // ones included (reordered back-to-front, but still drawn) — there is no toggle to turn
            // that off. Showing "a volume as its boundary" therefore means computing that boundary
            // ourselves (see volume_boundary_triangles()) and drawing it as a plain surface.
            auto *volume = polyscope::registerSurfaceMesh(GEOM_VOLUMES, vertices, volume_boundary);
            volume->setSurfaceColor(to_glm(m_config.geometry_volume_color));
            volume->setTransparency(static_cast<float>(m_config.model_transparency));
            volume->setEnabled(m_config.show_geometry_volumes);
        }
        if (!triangles.empty()) {
            auto *surface = polyscope::registerSurfaceMesh(GEOM_SURFACES, vertices, triangles);
            // See-through by default: a blocking is built around and inside its model, and an
            // opaque model hides the very thing being edited.
            surface->setSurfaceColor(to_glm(m_config.geometry_surface_color));
            surface->setTransparency(static_cast<float>(m_config.model_transparency));
            // Not necessarily redundant with "volumes" even when both are present: volumes shows
            // only the true outer skin of the tet mesh, while surfaces is every *tagged* surface —
            // on a model with an internal partition between two volumes (two_cubes.msh: 286 tagged
            // triangles against a 260-triangle outer boundary), that partition is real surface
            // content with no counterpart in "volumes" at all. So both follow their own config flag
            // independently; where they do coincide (a model with no internal walls, e.g. a single
            // solid), the two can be told apart in the Scene panel same as any other overlap.
            surface->setEnabled(m_config.show_geometry_surfaces);
        }
        if (!curve_segments.empty()) {
            auto *curves = polyscope::registerCurveNetwork(GEOM_CURVES, curve_points, curve_segments);
            curves->setColor(to_glm(m_config.geometry_curve_color));
            curves->setRadius(m_config.geometry_curve_radius);
            curves->setEnabled(m_config.show_geometry_curves);
        }
        if (!points.empty()) {
            auto *cloud = polyscope::registerPointCloud(GEOM_POINTS, points);
            cloud->setPointColor(to_glm(m_config.geometry_point_color));
            cloud->setPointRadius(m_config.geometry_point_radius);
            cloud->setEnabled(m_config.show_geometry_points);
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
        const auto vertices_shown = enabled_state(points_or_null(BLOCK_VERTICES));
        auto *cloud = polyscope::registerPointCloud(BLOCK_VERTICES, corners);
        cloud->setPointRadius(m_config.corner_radius);
        restore_enabled(cloud, vertices_shown, true);

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

        const auto edges_shown = enabled_state(curves_or_null(BLOCK_EDGES));
        const auto edge_points = m_blocking->edge_vertices(m_edge_samples);
        if (!edge_points.empty()) {
            auto *edges =
                polyscope::registerCurveNetwork(BLOCK_EDGES, edge_points, m_blocking->edge_segments(m_edge_samples));
            edges->setRadius(m_config.block_edge_radius);
            edges->setColor(to_glm(m_config.block_edge_color));
            restore_enabled(edges, edges_shown, m_show_block_edges);
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

        refresh_control_nets();

        // Block faces get their own structure rather than riding on the generated mesh: to_mesh()
        // only emits quads for standalone 2D blocks, so a hex block's 6 bounding faces would
        // otherwise be invisible — and unclassifiable to the eye. Sampled at the blocking's own
        // display resolution, not `m_subdivisions`, which now belongs to the Mesh section alone.
        const auto faces_shown = enabled_state(surface_or_null(BLOCK_FACES));
        const auto face_points = m_blocking->face_grid_vertices(m_edge_samples);
        if (!face_points.empty()) {
            auto *faces =
                polyscope::registerSurfaceMesh(BLOCK_FACES, face_points, m_blocking->face_grid_quads(m_edge_samples));
            const auto owners = m_blocking->face_grid_owners(m_edge_samples);
            apply_classification_colors(
                faces, m_blocking->face_classification_dims(), [&owners](auto *s, const auto &c) {
                    std::vector<glm::vec3> per_quad;
                    per_quad.reserve(owners.size());
                    for (const int owner : owners) {
                        per_quad.push_back(c[static_cast<std::size_t>(owner)]);
                    }
                    return s->addFaceColorQuantity("classification", per_quad);
                });
            restore_enabled(faces, faces_shown, true);
        } else if (polyscope::hasSurfaceMesh(BLOCK_FACES)) {
            polyscope::removeStructure(polyscope::getSurfaceMesh(BLOCK_FACES));
        }

        // The blocks themselves: one hex per 3-cell, which is exactly the mesh at subdivision 1.
        const auto blocks_shown = enabled_state(volume_or_null(BLOCK_BLOCKS));
        const auto blocks = m_blocking->mesh_hexes(1);
        if (!blocks.empty()) {
            restore_enabled(
                polyscope::registerHexMesh(BLOCK_BLOCKS, m_blocking->mesh_vertices(1), blocks), blocks_shown, true);
        } else if (polyscope::hasVolumeMesh(BLOCK_BLOCKS)) {
            polyscope::removeStructure(polyscope::getVolumeMesh(BLOCK_BLOCKS));
        }

        refresh_mesh();
    }

    void BiyApp::refresh_mesh() {
        // The mesh the blocking generates, at the chosen subdivision — a different object from the
        // block structure above, which is why `subdivisions` drives this and nothing else.
        const auto vertices = m_blocking->mesh_vertices(m_subdivisions);
        const auto quads = m_blocking->mesh_quads(m_subdivisions);
        const auto hexes = m_blocking->mesh_hexes(m_subdivisions);

        const auto quads_shown = enabled_state(surface_or_null(MESH_QUADS));
        if (!quads.empty()) {
            auto *mesh = polyscope::registerSurfaceMesh(MESH_QUADS, vertices, quads);
            mesh->setSurfaceColor(to_glm(m_config.mesh_color));
            restore_enabled(mesh, quads_shown, m_show_mesh);
        } else if (polyscope::hasSurfaceMesh(MESH_QUADS)) {
            polyscope::removeStructure(polyscope::getSurfaceMesh(MESH_QUADS));
        }

        const auto hexes_shown = enabled_state(volume_or_null(MESH_HEXES));
        if (!hexes.empty()) {
            auto *mesh = polyscope::registerHexMesh(MESH_HEXES, vertices, hexes);
            mesh->setColor(to_glm(m_config.mesh_color));
            restore_enabled(mesh, hexes_shown, m_show_mesh);
        } else if (polyscope::hasVolumeMesh(MESH_HEXES)) {
            polyscope::removeStructure(polyscope::getVolumeMesh(MESH_HEXES));
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

    void BiyApp::refresh_control_nets() {
        // At order 1 every control point coincides with a block corner (an edge's 2 are its
        // endpoints, a face's 4 its corners, a block's 8 its own), so there is nothing to reveal
        // and no structure is created.
        const bool curved = m_blocking && m_order > 1;
        const std::vector<std::array<double, 3>> none;
        const std::vector<std::array<int, 2>> no_segments;

        refresh_control_net(EDGE_CONTROL_POINTS,
                            EDGE_CONTROL_NET,
                            curved ? m_blocking->edge_control_points() : none,
                            curved ? m_blocking->edge_control_polygons() : no_segments,
                            m_config.edge_control_color,
                            m_show_edge_control);
        refresh_control_net(FACE_CONTROL_POINTS,
                            FACE_CONTROL_NET,
                            curved ? m_blocking->face_control_points() : none,
                            curved ? m_blocking->face_control_nets() : no_segments,
                            m_config.face_control_color,
                            m_show_face_control);
        refresh_control_net(BLOCK_CONTROL_POINTS,
                            BLOCK_CONTROL_NET,
                            curved ? m_blocking->block_control_points() : none,
                            curved ? m_blocking->block_control_lattices() : no_segments,
                            m_config.block_control_color,
                            m_show_block_control);
    }

    void BiyApp::refresh_control_net(const char *points_name,
                                     const char *net_name,
                                     const std::vector<std::array<double, 3>> &points,
                                     const std::vector<std::array<int, 2>> &segments,
                                     const std::array<float, 3> &color,
                                     bool visible) {
        if (points.empty()) {
            if (polyscope::hasPointCloud(points_name)) {
                polyscope::removeStructure(polyscope::getPointCloud(points_name));
            }
            if (polyscope::hasCurveNetwork(net_name)) {
                polyscope::removeStructure(polyscope::getCurveNetwork(net_name));
            }
            return;
        }

        // Points and segments are 2 structures rather than a single curve network (which would draw
        // both) so the points can be spheres big enough to catch the eye while the net stays a thin
        // scaffold behind them.
        auto *cloud = polyscope::registerPointCloud(points_name, points);
        cloud->setPointRadius(m_config.control_point_radius);
        cloud->setPointColor(to_glm(color));
        cloud->setEnabled(visible);

        auto *net = polyscope::registerCurveNetwork(net_name, points, segments);
        net->setRadius(m_config.control_polygon_radius);
        net->setColor(to_glm(color));
        net->setEnabled(visible);
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

        draw_operations_panel();
        draw_scene_panel();
        // Kept from Polyscope's own set: a separate window reporting what was last clicked, which
        // is also how corners get picked up. Its main panel and structure list are replaced above.
        polyscope::buildPickGui();

        handle_drag();
    }

    void BiyApp::draw_operations_panel() {
        ImGui::SetNextWindowPos(ImVec2(polyscope::internal::imguiStackMargin, polyscope::internal::imguiStackMargin));
        ImGui::SetNextWindowSize(ImVec2(polyscope::internal::leftWindowsWidth, 0.0f));
        ImGui::Begin("BIY operations", nullptr);

        // The parts of Polyscope's own panel worth keeping, rebuilt here from its public builders
        // rather than left in a window of their own. Its "Debug" section is dropped: it is an
        // internal texture inspector, inlined in buildPolyscopeGui() and not separately callable.
        if (ImGui::Button("Reset view")) {
            polyscope::view::flyToHomeView();
        }
        ImGui::SameLine();
        if (ImGui::Button("Screenshot")) {
            polyscope::ScreenshotOptions options;
            options.transparentBackground = polyscope::options::screenshotTransparency;
            polyscope::screenshot(options);
        }
        polyscope::view::buildViewGui();
        polyscope::render::engine->buildEngineGui();
        ImGui::Separator();

        draw_panel();
        ImGui::End();
    }

    void BiyApp::draw_scene_entry(const char *label, polyscope::Structure *structure) {
        if (structure == nullptr) {
            ImGui::BeginDisabled();
            ImGui::BulletText("%s (none)", label);
            ImGui::EndDisabled();
            return;
        }
        // Delegated to the structure itself, so the visibility controls are Polyscope's real ones
        // rather than a reimplementation of them.
        ImGui::PushID(label);
        structure->buildUI();
        ImGui::PopID();
    }

    void BiyApp::draw_scene_panel() {
        // On the opposite side from "BIY operations" — its own column, not stacked under the other
        // one — so the two read as independent panels rather than one long list.
        //
        // A fixed strip at the bottom is reserved for Polyscope's own "Selection" box
        // (polyscope::buildPickGui(), drawn right after this), whether or not something is
        // currently picked: sizing this panel around the *presence* of a selection would make it
        // grow and shrink as picks come and go, which reads as flicker more than as useful space.
        // buildPickGui() finds where to start via internal::lastRightSideFreeY, which normally only
        // Polyscope's own stacked panels update — since biy draws no such panel above it, that stays
        // near 0 and Selection would otherwise land right on top of this one.
        constexpr float SELECTION_RESERVED_HEIGHT = 160.0f;
        const float margin = polyscope::internal::imguiStackMargin;
        const float width = polyscope::internal::rightWindowsWidth;
        const float height =
            polyscope::view::windowHeight - 3 * margin - SELECTION_RESERVED_HEIGHT * polyscope::options::uiScale;

        ImGui::SetNextWindowPos(ImVec2(polyscope::view::windowWidth - width - margin, margin));
        ImGui::SetNextWindowSize(ImVec2(width, height));
        ImGui::Begin("Scene", nullptr);

        if (ImGui::CollapsingHeader("Geometry", ImGuiTreeNodeFlags_DefaultOpen)) {
            draw_scene_entry("volumes", surface_or_null(GEOM_VOLUMES));
            draw_scene_entry("surfaces", surface_or_null(GEOM_SURFACES));
            draw_scene_entry("curves", curves_or_null(GEOM_CURVES));
            draw_scene_entry("points", points_or_null(GEOM_POINTS));
        }

        if (ImGui::CollapsingHeader("Blocking", ImGuiTreeNodeFlags_DefaultOpen)) {
            draw_scene_entry("blocks", volume_or_null(BLOCK_BLOCKS));
            draw_scene_entry("faces", surface_or_null(BLOCK_FACES));
            draw_scene_entry("edges", curves_or_null(BLOCK_EDGES));
            draw_scene_entry("vertices", points_or_null(BLOCK_VERTICES));
        }

        if (ImGui::CollapsingHeader("Mesh")) {
            draw_scene_entry("quads", surface_or_null(MESH_QUADS));
            draw_scene_entry("hexes", volume_or_null(MESH_HEXES));
        }

        ImGui::End();

        // See this panel's own reserved-strip comment above: this is what actually makes
        // buildPickGui() start below the strip rather than back at the top of the column.
        polyscope::internal::lastRightSideFreeY = margin + height;
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
        // The step sizes are given explicitly because InputInt's own default fast step is 100 —
        // meaningless over a 1..20 range, and reachable without meaning to: ImGui applies it when
        // io.KeyCtrl is set, which on macOS means Cmd (ConfigMacOSXBehaviors swaps the two). With
        // the cap in place that turned every arrow click into a jump straight to one end or the
        // other.
        if (ImGui::InputInt("subdivisions", &m_subdivisions, 1, 5)) {
            m_subdivisions = std::clamp(m_subdivisions, 1, MAX_SUBDIVISIONS);
            refresh_view();
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("1 to %d (Cmd-click the arrows to step by 5). Each step costs cubically on a "
                              "3D blocking — %d subdivisions already means %d hexes per block.",
                              MAX_SUBDIVISIONS,
                              MAX_SUBDIVISIONS,
                              MAX_SUBDIVISIONS * MAX_SUBDIVISIONS * MAX_SUBDIVISIONS);
        }

        if (ImGui::Checkbox("Show block edges", &m_show_block_edges)) refresh_view();

        ImGui::BeginDisabled(m_order <= 1);
        ImGui::TextUnformatted("Control points");
        if (ImGui::Checkbox("edges", &m_show_edge_control)) refresh_view();
        ImGui::SameLine();
        if (ImGui::Checkbox("faces", &m_show_face_control)) refresh_view();
        ImGui::SameLine();
        if (ImGui::Checkbox("blocks", &m_show_block_control)) refresh_view();
        ImGui::EndDisabled();
        if (m_order <= 1 && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
            ImGui::SetTooltip("At order 1 every control point is just a block corner.");
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
            if (pick.isHit && pick.structureName == BLOCK_VERTICES) {
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
