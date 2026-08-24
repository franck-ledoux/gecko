#include "BiyApp.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

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
        /** @brief Cut-mode preview: the sheet under the cursor, and where the cut would land on each
         * of its edges. Separate short-lived structures, registered only while hovering. */
        constexpr const char *SHEET_EDGES = "sheet";
        constexpr const char *CUT_POINTS = "cut";
        /** @brief Delete-mode preview: the one block a click would remove. */
        constexpr const char *DELETE_PREVIEW = "to delete";

        glm::vec3 to_glm(const std::array<float, 3> &c) { return {c[0], c[1], c[2]}; }

        /**
         * @brief The worst departure from straightness anywhere in the blocking, as a line for the
         * terminal.
         *
         * Reported after every cut and deletion, because a straight blocking has to stay straight
         * through both — cutting a straight edge gives straight halves and deleting moves nothing —
         * and because when it does not, this number says *where the fault is*: a large value means
         * the geometry itself was corrupted, a value at rounding while the screen shows curves means
         * the drawing was, and those are opposite places to look.
         */
        std::string bend_report(app::BlockingFacade &blocking) {
            const auto bends = blocking.edge_bends();
            if (bends.empty()) return {};
            const auto worst = std::max_element(bends.begin(), bends.end());
            std::ostringstream out;
            out << std::scientific << std::setprecision(2) << " worst edge bend " << *worst << " (edge "
                << std::distance(bends.begin(), worst) << ")";
            return out.str();
        }

        /** @brief A cut parameter, short enough to sit in the status line. */
        std::string format_param(double value) {
            std::ostringstream out;
            out << std::fixed << std::setprecision(3) << value;
            return out.str();
        }

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
        : m_model(std::make_unique<app::GeomModelFacade>(model_path)), m_model_path(model_path), m_order(order) {
        std::string config_message;
        m_config = BiyConfig::load("biy_config.json", config_message);
        m_show_block_edges = m_config.show_block_edges;
        m_show_mesh = m_config.show_mesh;
        m_mesh_coloring = m_config.mesh_color_by_block ? MeshColoring::ByBlock : MeshColoring::Uniform;
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
        if (mode != MouseMode::Cut && mode != MouseMode::Collapse) {
            m_hover_edge.reset();
            m_sheet.clear();
            m_last_cut_mouse = glm::vec2(-1.0f, -1.0f);
            refresh_cut_preview();
        }
        if (mode != MouseMode::Delete) {
            m_hover_block.reset();
            m_last_delete_mouse = glm::vec2(-1.0f, -1.0f);
            refresh_delete_preview();
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

    app::BlockingFacade &BiyApp::blocking() {
        if (!m_blocking) m_blocking = std::make_unique<app::BlockingFacade>(*m_model, m_order);
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

    int BiyApp::display_position(int ADim, int AId) const {
        const auto ids = (ADim == 1) ? m_blocking->edge_ids() : m_blocking->block_ids();
        const auto found = std::find(ids.begin(), ids.end(), AId);
        return (found == ids.end()) ? -1 : static_cast<int>(std::distance(ids.begin(), found));
    }

    glm::vec3 BiyApp::block_color(int index) {
        // Hues stepped by the golden ratio rather than picked from a fixed palette: a blocking has
        // no bound on how many blocks it holds, and a palette of N runs out and starts repeating —
        // 2 neighbours sharing a colour reads as one block, which is the one thing this must not do.
        // Stepping by an irrational fraction of the circle keeps consecutive indices far apart and
        // never lands twice on the same hue, however many there are.
        constexpr float golden = 0.6180339887498949f;
        const float hue = std::fmod(static_cast<float>(index) * golden, 1.0f);
        // Held well below full saturation and brightness so the mesh still reads as a mesh under
        // Polyscope's shading, and so its edges stay visible against it.
        constexpr float saturation = 0.55f;
        constexpr float value = 0.92f;

        const float h = hue * 6.0f;
        const int sector = static_cast<int>(h) % 6;
        const float f = h - std::floor(h);
        const float p = value * (1.0f - saturation);
        const float q = value * (1.0f - saturation * f);
        const float t = value * (1.0f - saturation * (1.0f - f));
        switch (sector) {
            case 0:
                return {value, t, p};
            case 1:
                return {q, value, p};
            case 2:
                return {p, value, t};
            case 3:
                return {p, q, value};
            case 4:
                return {t, p, value};
            default:
                return {value, p, q};
        }
    }

    void BiyApp::refresh_mesh() {
        // The mesh the blocking generates, at the chosen subdivision — a different object from the
        // block structure above, which is why `subdivisions` drives this and nothing else.
        const auto vertices = m_blocking->mesh_vertices(m_subdivisions);
        const auto quads = m_blocking->mesh_quads(m_subdivisions);
        const auto hexes = m_blocking->mesh_hexes(m_subdivisions);

        // One colour per block, spread over every cell that block generated. Built here rather
        // than left to a Polyscope colormap: the owner indices are block *identities*, not a scale,
        // and a colormap would shade them as if 0 and 1 were closer than 0 and 9.
        const bool by_block = m_mesh_coloring == MeshColoring::ByBlock;
        const auto per_cell = [this](const std::vector<int> &owners) {
            std::vector<glm::vec3> colors;
            colors.reserve(owners.size());
            for (const int owner : owners) {
                colors.push_back(block_color(owner));
            }
            return colors;
        };

        const auto quads_shown = enabled_state(surface_or_null(MESH_QUADS));
        if (!quads.empty()) {
            auto *mesh = polyscope::registerSurfaceMesh(MESH_QUADS, vertices, quads);
            mesh->setSurfaceColor(to_glm(m_config.mesh_color));
            if (by_block) {
                mesh->addFaceColorQuantity("block", per_cell(m_blocking->mesh_quad_owners(m_subdivisions)))
                    ->setEnabled(true);
            }
            restore_enabled(mesh, quads_shown, m_show_mesh);
        } else if (polyscope::hasSurfaceMesh(MESH_QUADS)) {
            polyscope::removeStructure(polyscope::getSurfaceMesh(MESH_QUADS));
        }

        const auto hexes_shown = enabled_state(volume_or_null(MESH_HEXES));
        if (!hexes.empty()) {
            auto *mesh = polyscope::registerHexMesh(MESH_HEXES, vertices, hexes);
            mesh->setColor(to_glm(m_config.mesh_color));
            if (by_block) {
                mesh->addCellColorQuantity("block", per_cell(m_blocking->mesh_hex_owners(m_subdivisions)))
                    ->setEnabled(true);
            }
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
        draw_gizmo();

        handle_drag();
        handle_cut();
        handle_collapse();
        handle_delete();
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

    void BiyApp::draw_gizmo() {
        // A flat 2D compass rather than a real 3D widget, deliberately: "bottom-right of the
        // screen" is a screen-space request, and a world-space object anchored to the camera would
        // need careful per-frame offset/scale math to stay pinned there under every zoom level
        // without ever clipping through nearby scene geometry. Projecting the 6 world axes through
        // the camera's own right/up/look basis (getCameraFrame()) onto a small ImGui canvas gets the
        // same "always shows current orientation" property with none of that risk, at the cost of
        // not being real depth-tested geometry — approximated below by depth-sorting and shrinking
        // the farther dots instead.
        const float scale = polyscope::options::uiScale;
        const float size = static_cast<float>(m_config.gizmo_size) * scale;
        const float margin = polyscope::internal::imguiStackMargin;
        const ImVec2 window_pos(polyscope::view::windowWidth - size - margin,
                                polyscope::view::windowHeight - size - margin);

        ImGui::SetNextWindowPos(window_pos);
        ImGui::SetNextWindowSize(ImVec2(size, size));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.0f, 0.0f, 0.0f, 0.0f));
        ImGui::Begin("##orientation_gizmo",
                     nullptr,
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoScrollbar |
                         ImGuiWindowFlags_NoScrollWithMouse | ImGuiWindowFlags_NoCollapse |
                         ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing |
                         ImGuiWindowFlags_NoNav | ImGuiWindowFlags_NoBackground);

        const auto [look, up, right] = polyscope::view::getCameraParametersForCurrentView().getCameraFrame();
        const ImVec2 center(window_pos.x + size * 0.5f, window_pos.y + size * 0.5f);
        const float radius = size * 0.35f;
        const float dot_radius = static_cast<float>(m_config.gizmo_dot_radius) * scale;

        struct Axis {
            glm::vec3 dir;
            const char *label;
            const std::array<float, 3> &color;
            bool positive;
        };
        const std::array<Axis, 6> axes{{{{1, 0, 0}, "X", m_config.gizmo_color_x, true},
                                        {{-1, 0, 0}, "", m_config.gizmo_color_x, false},
                                        {{0, 1, 0}, "Y", m_config.gizmo_color_y, true},
                                        {{0, -1, 0}, "", m_config.gizmo_color_y, false},
                                        {{0, 0, 1}, "Z", m_config.gizmo_color_z, true},
                                        {{0, 0, -1}, "", m_config.gizmo_color_z, false}}};

        // Farthest-from-camera first, so nearer dots paint over farther ones where they overlap —
        // the same painter's-algorithm depth cue a real 3D gizmo gets from the depth buffer.
        std::array<int, 6> order{0, 1, 2, 3, 4, 5};
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return glm::dot(axes[a].dir, look) > glm::dot(axes[b].dir, look);
        });

        ImDrawList *draw_list = ImGui::GetWindowDrawList();
        for (const int i : order) {
            const Axis &axis = axes[i];
            const float depth = glm::dot(axis.dir, look); // >0 into the screen, <0 toward the viewer
            const ImVec2 dot_pos(center.x + glm::dot(axis.dir, right) * radius,
                                 center.y - glm::dot(axis.dir, up) * radius);
            const float dot_size = dot_radius * (1.15f - 0.35f * depth);
            const ImU32 color =
                ImGui::ColorConvertFloat4ToU32(ImVec4(axis.color[0], axis.color[1], axis.color[2], 1.0f));

            ImGui::SetCursorScreenPos(ImVec2(dot_pos.x - dot_size, dot_pos.y - dot_size));
            ImGui::PushID(i);
            ImGui::InvisibleButton("##axis", ImVec2(2 * dot_size, 2 * dot_size));
            if (ImGui::IsItemClicked()) {
                fly_to_axis_view(axis.dir);
            }
            ImGui::PopID();

            if (axis.positive) {
                draw_list->AddCircleFilled(dot_pos, dot_size, color);
                const ImVec2 label_size = ImGui::CalcTextSize(axis.label);
                draw_list->AddText(ImVec2(dot_pos.x - label_size.x * 0.5f, dot_pos.y - label_size.y * 0.5f),
                                   IM_COL32_BLACK,
                                   axis.label);
            } else {
                draw_list->AddCircle(dot_pos, dot_size, color, 0, 1.5f * scale);
            }
        }

        ImGui::End();
        ImGui::PopStyleColor();
    }

    void BiyApp::fly_to_axis_view(const glm::vec3 &axis_dir) {
        const glm::vec3 target = polyscope::state::center();
        const float distance = 1.5f * polyscope::state::lengthScale;
        // Same degenerate case computeHomeView() itself guards against: lookAt() needs an up vector
        // that isn't parallel to the view direction. Looking along the scene's own up axis is
        // exactly that case (e.g. the +Y view when Y is up), so the scene's front axis stands in for
        // up there instead.
        const glm::vec3 scene_up = polyscope::view::getUpVec();
        const bool degenerate = std::abs(glm::dot(glm::normalize(axis_dir), glm::normalize(scene_up))) > 0.99f;
        const glm::vec3 up_for_view = degenerate ? polyscope::view::getFrontVec() : scene_up;

        polyscope::view::lookAt(target + axis_dir * distance, target, up_for_view, true);
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
            if (ImGui::IsKeyPressed(ImGuiKey_X)) set_mouse_mode(MouseMode::Cut);
            if (ImGui::IsKeyPressed(ImGuiKey_S)) set_mouse_mode(MouseMode::Collapse);
            if (ImGui::IsKeyPressed(ImGuiKey_D)) set_mouse_mode(MouseMode::Delete);
        }

        ImGui::TextUnformatted("Mouse mode");
        if (ImGui::RadioButton("Camera (C)", m_mode == MouseMode::Camera)) set_mouse_mode(MouseMode::Camera);
        ImGui::SameLine();
        if (ImGui::RadioButton("Edit (E)", m_mode == MouseMode::Edit)) set_mouse_mode(MouseMode::Edit);
        ImGui::SameLine();
        if (ImGui::RadioButton("Cut (X)", m_mode == MouseMode::Cut)) set_mouse_mode(MouseMode::Cut);
        ImGui::SameLine();
        if (ImGui::RadioButton("Collapse (S)", m_mode == MouseMode::Collapse)) set_mouse_mode(MouseMode::Collapse);
        ImGui::SameLine();
        if (ImGui::RadioButton("Delete (D)", m_mode == MouseMode::Delete)) set_mouse_mode(MouseMode::Delete);
        switch (m_mode) {
            case MouseMode::Edit:
                ImGui::TextWrapped("Drag a block corner to move it. Camera navigation is off.");
                break;
            case MouseMode::Cut:
                ImGui::TextWrapped("Point at a block edge: the whole sheet that would be cut lights up, and the "
                                   "markers show where. Click, or press Space, to cut it. Camera navigation "
                                   "is off.");
                break;
            case MouseMode::Collapse:
                ImGui::TextWrapped("Point at a block edge: the whole sheet lights up. Click, or press Space, to "
                                   "take that whole layer out and glue back what was either side of it. Where 2 "
                                   "corners merge, the one on the more constrained bit of the model wins — and "
                                   "it is refused outright when that would merge 2 different model vertices. "
                                   "Taking the last layer out empties the blocking. Camera navigation is off.");
                break;
            case MouseMode::Delete:
                ImGui::TextWrapped("Point at a block: it lights up. Click, or press Space, to delete it — along "
                                   "with everything that existed only because of it. Camera navigation is off.");
                break;
            default:
                ImGui::TextWrapped("Rotate/pan/zoom the view. Switch to Edit to move corners, Cut to split "
                                   "blocks, Collapse to take a whole layer out, Delete to remove one block.");
                break;
        }

        if (m_mode == MouseMode::Collapse) {
            if (m_hover_edge && !m_sheet.empty()) {
                ImGui::Text("Sheet: %zu edges", m_sheet.size());
            } else if (m_hover_edge) {
                ImGui::TextWrapped("This sheet closes back onto itself — it cannot be collapsed.");
            } else {
                ImGui::TextUnformatted("Sheet: none under the cursor");
            }
        }
        if (m_mode == MouseMode::Cut) {
            // The cursor drives this every time it moves; typing in it is the way to land a cut on
            // an exact value, which pointing cannot do.
            ImGui::SetNextItemWidth(120.0f * polyscope::options::uiScale);
            ImGui::SliderFloat("Cut at", &m_cut_param, 0.01f, 0.99f, "%.3f");
            if (ImGui::IsItemEdited() && m_hover_edge) refresh_cut_preview();
            if (m_hover_edge && !m_sheet.empty()) {
                ImGui::Text("Sheet: %zu edges", m_sheet.size());
            } else if (m_hover_edge) {
                ImGui::TextWrapped("This sheet closes back onto itself — no even cut exists.");
            } else {
                ImGui::TextUnformatted("Sheet: none under the cursor");
            }
        }
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

        ImGui::TextUnformatted("Mesh color");
        ImGui::SameLine();
        if (ImGui::RadioButton("uniform", m_mesh_coloring == MeshColoring::Uniform)) {
            m_mesh_coloring = MeshColoring::Uniform;
            refresh_mesh();
        }
        ImGui::SameLine();
        if (ImGui::RadioButton("per block", m_mesh_coloring == MeshColoring::ByBlock)) {
            m_mesh_coloring = MeshColoring::ByBlock;
            refresh_mesh();
        }

        ImGui::SetNextItemWidth(120.0f);
        if (ImGui::InputInt("order", &m_order, 1, 1)) {
            m_order = std::clamp(m_order, 1, MAX_ORDER);
            if (m_blocking) {
                // Not just more control points: the whole structure is refitted at the new order, so
                // an edge lying on a curved model curve stops being its chord and starts following
                // it. Uses the panel's own tolerances, the same ones Classify uses.
                m_blocking->set_degree(m_order, m_tol_vertex, m_tol_curve, m_tol_surface);
                refresh_view();
                m_status = "Rebuilt at order " + std::to_string(m_order) + bend_report(*m_blocking);
                std::cout << "biy [order " << m_order << "]: blocking refitted." << std::endl;
            }
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("1 to %d. 1 is straight edges. Changing it refits every edge, face and block onto "
                              "the model at the new order — topology and classification are kept.",
                              MAX_ORDER);
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

    void BiyApp::update_cut_hover(glm::vec2 screen_coords) {
        const auto clear = [this] {
            m_hover_edge.reset();
            m_sheet.clear();
        };

        const polyscope::PickResult pick = polyscope::pickAtScreenCoords(screen_coords);
        if (!pick.isHit) {
            clear();
            return;
        }

        // Aiming a cut means aiming at any of 3 things, not just the bare edge: the edge itself, the
        // highlight drawn over it, or one of the markers showing where the cut would land. All 3 sit
        // on the same curve, so all 3 resolve below to the same edge and to much the same point
        // along it — and accepting all 3 is what makes the target the width of a marker rather than
        // the width of a hairline.
        //
        // The other 2 are not a nicety. The highlight is drawn thicker than the edges it traces and
        // covers them, and each marker is thicker still and lands exactly under the cursor that
        // summoned it — so a moment after a sheet lights up the cursor is no longer picking the edge
        // at all. Reading that as "nothing there" made the preview flicker on and off as the cursor
        // moved along an edge, and a cut aimed through it do nothing.
        //
        // They go through the same search rather than freezing the hover, deliberately: a marker
        // follows the cursor, so holding the parameter while over one would pin the cut where it
        // first appeared and leave it unmovable.
        if (pick.structureName != BLOCK_EDGES && pick.structureName != SHEET_EDGES &&
            pick.structureName != CUT_POINTS) {
            clear();
            return;
        }

        // Polyscope reports where in space the click landed, which is all that is needed: the block
        // edges are drawn as polylines sampled from the very same curves the facade indexes, so
        // finding the nearest sampled segment recovers both which edge it was and how far along —
        // without depending on how Polyscope happens to number a curve network's elements.
        const auto points = m_blocking->edge_vertices(m_edge_samples);
        const auto per_edge = static_cast<std::size_t>(m_edge_samples) + 1;
        if (points.empty() || points.size() % per_edge != 0) {
            clear();
            return;
        }

        const glm::vec3 target = pick.position;
        double best = std::numeric_limits<double>::max();
        int best_edge = -1;
        double best_param = 0.5;
        for (std::size_t e = 0; e * per_edge < points.size(); ++e) {
            for (std::size_t k = 0; k + 1 < per_edge; ++k) {
                const auto &a = points[e * per_edge + k];
                const auto &b = points[e * per_edge + k + 1];
                const glm::vec3 pa(static_cast<float>(a[0]), static_cast<float>(a[1]), static_cast<float>(a[2]));
                const glm::vec3 pb(static_cast<float>(b[0]), static_cast<float>(b[1]), static_cast<float>(b[2]));
                const glm::vec3 ab = pb - pa;
                const float len2 = glm::dot(ab, ab);
                const float f = (len2 > 0.0f) ? glm::clamp(glm::dot(target - pa, ab) / len2, 0.0f, 1.0f) : 0.0f;
                const double distance = glm::length(target - (pa + ab * f));
                if (distance < best) {
                    best = distance;
                    best_edge = static_cast<int>(e);
                    best_param = (static_cast<double>(k) + f) / static_cast<double>(m_edge_samples);
                }
            }
        }
        if (best_edge < 0) {
            clear();
            return;
        }

        // Cutting in half is far and away the common case, so the middle gets a magnet rather than
        // asking for pixel-perfect aiming.
        if (std::abs(best_param - 0.5) < m_config.cut_snap_tolerance) best_param = 0.5;

        // The search ran over the display arrays, so `best_edge` is a position; what is kept is
        // the id at that position, since the hover outlives the frame it was made in. The sheet is
        // the other way round: it is drawn, so it stays positional.
        m_hover_edge = m_blocking->edge_ids()[static_cast<std::size_t>(best_edge)];
        m_cut_param = static_cast<float>(best_param);
        m_sheet = m_blocking->sheet_edges(*m_hover_edge);
    }

    void BiyApp::refresh_cut_preview() {
        const bool cutting = m_mode == MouseMode::Cut;
        const bool show = (cutting || m_mode == MouseMode::Collapse) && m_hover_edge.has_value() && !m_sheet.empty();
        if (!show) {
            if (polyscope::hasCurveNetwork(SHEET_EDGES)) {
                polyscope::removeStructure(polyscope::getCurveNetwork(SHEET_EDGES));
            }
            if (polyscope::hasPointCloud(CUT_POINTS)) {
                polyscope::removeStructure(polyscope::getPointCloud(CUT_POINTS));
            }
            return;
        }

        // Rebuilt as its own compact point list rather than reusing the full edge_vertices() array:
        // Polyscope draws a sphere per curve-network node, so carrying along every point of every
        // unhighlighted edge would bury the sheet under thousands of them.
        const auto points = m_blocking->edge_vertices(m_edge_samples);
        const auto per_edge = static_cast<std::size_t>(m_edge_samples) + 1;
        std::vector<std::array<double, 3>> sheet_points;
        std::vector<std::array<int, 2>> sheet_segments;
        for (const int e : m_sheet) {
            const auto base = static_cast<int>(sheet_points.size());
            for (std::size_t k = 0; k < per_edge; ++k) {
                sheet_points.push_back(points[static_cast<std::size_t>(e) * per_edge + k]);
            }
            for (int k = 0; k < m_edge_samples; ++k) {
                sheet_segments.push_back({base + k, base + k + 1});
            }
        }

        const auto &color = cutting ? m_config.sheet_color : m_config.collapse_color;
        auto *net = polyscope::registerCurveNetwork(SHEET_EDGES, sheet_points, sheet_segments);
        net->setColor(glm::vec3(color[0], color[1], color[2]));
        net->setRadius(static_cast<float>(m_config.sheet_radius));

        // Collapsing has no parameter, so it has nothing to mark: the markers say "the cut lands
        // here", and there is no here.
        if (!cutting) {
            if (polyscope::hasPointCloud(CUT_POINTS)) {
                polyscope::removeStructure(polyscope::getPointCloud(CUT_POINTS));
            }
            return;
        }

        const auto cut_points = m_blocking->sheet_cut_points(*m_hover_edge, m_cut_param);
        auto *cloud = polyscope::registerPointCloud(CUT_POINTS, cut_points);
        cloud->setPointColor(
            glm::vec3(m_config.cut_point_color[0], m_config.cut_point_color[1], m_config.cut_point_color[2]));
        cloud->setPointRadius(static_cast<float>(m_config.cut_point_radius));
    }

    void BiyApp::handle_cut() {
        if (!m_blocking || m_mode != MouseMode::Cut) return;
        ImGuiIO &io = ImGui::GetIO();

        // Space cuts too. A trackpad tap is a genuinely unreliable way to say "here": it fires a
        // click and a small cursor jolt together, and it can be shorter than a frame. A key is
        // neither, and it also lets a careful cut be placed without the pointer being nudged off the
        // edge in the act of pressing.
        const bool by_key = !io.WantCaptureKeyboard && ImGui::IsKeyPressed(ImGuiKey_Space);
        const bool by_click = !io.WantCaptureMouse && io.MouseClicked[0];

        if (by_click || by_key) {
            // Acts on the sheet the preview is *currently showing*, without re-testing what is under
            // the cursor first. That is the point: the pointer often shifts a pixel or two in the act
            // of clicking, and on edges this thin that is enough to fall off them — re-testing here
            // would clear the hover and swallow the click, which is exactly the "sometimes nothing
            // happens" this fixes. What you saw highlighted is what gets cut.
            perform_cut(by_key ? "space" : "click");
            return;
        }
        if (io.WantCaptureMouse) return;

        const glm::vec2 mouse{io.MousePos.x, io.MousePos.y};
        // Only re-tested when the cursor actually moved: a hover test is a Polyscope pick, and a
        // pick is a render pass — running one every frame for a still cursor would cost real
        // framerate for no new information.
        if (mouse != m_last_cut_mouse) {
            m_last_cut_mouse = mouse;
            update_cut_hover(mouse);
            refresh_cut_preview();
        }
    }

    void BiyApp::perform_cut(const char *trigger) {
        // Everything here also goes to the terminal the Python console runs in: the status line is
        // one line at the bottom of a panel and is easy to miss, and "I pressed it and nothing
        // happened" needs an answer that stays on screen.
        const auto report = [trigger](const std::string &message) {
            std::cout << "biy [cut by " << trigger << "]: " << message << std::endl;
        };

        if (!m_hover_edge) {
            m_status = "Nothing under the cursor to cut";
            report("no block edge under the cursor. Point at one — the sheet lights up when you are on "
                   "it — then click or press space. If the edges are hidden in the Scene panel, "
                   "nothing can be picked.");
            return;
        }
        if (m_sheet.empty()) {
            m_status = "That sheet cannot be cut evenly";
            report("that sheet closes back onto itself, so no single cut splits it evenly. Refused "
                   "rather than cut somewhere arbitrary.");
            return;
        }

        const std::size_t cut_edges = m_sheet.size();
        const std::size_t before = m_blocking->nb_cells(3);
        if (!m_blocking->cut_sheet(*m_hover_edge, m_cut_param)) {
            m_status = "Cut refused at t=" + format_param(m_cut_param);
            report("the kernel refused the cut at t=" + format_param(m_cut_param) + " on edge " +
                   std::to_string(*m_hover_edge) + ".");
            return;
        }
        const std::size_t after = m_blocking->nb_cells(3);

        // The whole preview describes edges that no longer exist, and the indices behind it have all
        // shifted, so it goes before anything is drawn again.
        m_hover_edge.reset();
        m_sheet.clear();
        m_last_cut_mouse = glm::vec2(-1.0f, -1.0f);
        refresh_cut_preview();
        refresh_view();

        m_status = "Cut " + std::to_string(cut_edges) + " edges at t=" + format_param(m_cut_param) +
                   (after > before ? " — blocks " + std::to_string(before) + " to " + std::to_string(after) : "");
        report("cut " + std::to_string(cut_edges) + " edges at t=" + format_param(m_cut_param) + ", blocks " +
               std::to_string(before) + " -> " + std::to_string(after) + "." + bend_report(*m_blocking));
    }

    void BiyApp::handle_collapse() {
        if (!m_blocking || m_mode != MouseMode::Collapse) return;
        ImGuiIO &io = ImGui::GetIO();

        // Space as well as a click, and the sheet under the cursor read only when the cursor moves —
        // both for the reasons `handle_cut()` spells out.
        const bool by_key = !io.WantCaptureKeyboard && ImGui::IsKeyPressed(ImGuiKey_Space);
        const bool by_click = !io.WantCaptureMouse && io.MouseClicked[0];

        if (by_click || by_key) {
            perform_collapse(by_key ? "space" : "click");
            return;
        }
        if (io.WantCaptureMouse) return;

        const glm::vec2 mouse{io.MousePos.x, io.MousePos.y};
        if (mouse != m_last_cut_mouse) {
            m_last_cut_mouse = mouse;
            update_cut_hover(mouse);
            refresh_cut_preview();
        }
    }

    void BiyApp::perform_collapse(const char *trigger) {
        const auto report = [trigger](const std::string &message) {
            std::cout << "biy [collapse by " << trigger << "]: " << message << std::endl;
        };

        if (!m_hover_edge) {
            m_status = "Nothing under the cursor to collapse";
            report("no block edge under the cursor. Point at one — the sheet lights up when you are "
                   "on it — then click or press space.");
            return;
        }
        if (m_sheet.empty()) {
            m_status = "That sheet cannot be collapsed";
            report("that sheet closes back onto itself, so it has no coherent 2 sides to glue "
                   "together. Refused rather than collapsed into something arbitrary.");
            return;
        }

        const std::size_t sheet_edges = m_sheet.size();
        const std::size_t before = m_blocking->nb_cells(3);
        if (!m_blocking->delete_sheet(*m_hover_edge, m_tol_vertex, m_tol_curve, m_tol_surface)) {
            m_status = "Collapse refused";
            report("the kernel refused to collapse the sheet through edge " + std::to_string(*m_hover_edge) +
                   ". On a classified blocking that is almost always one edge of the sheet joining 2 corners "
                   "that sit on 2 different vertices of the model: merging them would leave one of those "
                   "vertices with no corner of the block structure on it, and nothing else records where it "
                   "was. Otherwise it is a sheet closing back onto itself, or one that would leave an edge "
                   "joining a corner to itself.");
            return;
        }
        const std::size_t after = m_blocking->nb_cells(3);

        m_hover_edge.reset();
        m_sheet.clear();
        m_last_cut_mouse = glm::vec2(-1.0f, -1.0f);
        refresh_cut_preview();
        refresh_view();

        m_status = "Collapsed a sheet of " + std::to_string(sheet_edges) + " edges — blocks " + std::to_string(before) +
                   " to " + std::to_string(after);
        report("collapsed a sheet of " + std::to_string(sheet_edges) + " edges, blocks " + std::to_string(before) +
               " -> " + std::to_string(after) + "." + bend_report(*m_blocking));
    }

    void BiyApp::update_delete_hover(glm::vec2 screen_coords) {
        const auto clear = [this] { m_hover_block.reset(); };

        const polyscope::PickResult pick = polyscope::pickAtScreenCoords(screen_coords);
        if (!pick.isHit) {
            clear();
            return;
        }
        // The highlight is drawn over the block it marks and so is what the cursor picks from the
        // moment it appears. Reading that as "nothing there" would drop the aim the instant it was
        // made — the same trap the cut preview fell into. Here there is nothing to re-derive from
        // it, so holding the current block is exactly right.
        if (pick.structureName == DELETE_PREVIEW) return;
        if (pick.structureName != BLOCK_BLOCKS || !polyscope::hasVolumeMesh(BLOCK_BLOCKS)) {
            clear();
            return;
        }

        // The blocks are registered as one hex per 3-cell, in the blocking's own block order, so a
        // picked cell index *is* the index delete_block() speaks.
        const auto hit = polyscope::getVolumeMesh(BLOCK_BLOCKS)->interpretPickResult(pick);
        if (hit.elementType != polyscope::VolumeMeshElement::CELL) {
            clear();
            return;
        }
        // The blocks are drawn one hex per 3-cell in the blocking's own order, so a picked cell
        // index is a position; kept as the id at that position, for the reason `m_hover_edge` gives.
        m_hover_block = m_blocking->block_ids()[hit.index];
    }

    void BiyApp::refresh_delete_preview() {
        const bool show = m_mode == MouseMode::Delete && m_hover_block.has_value();
        if (!show) {
            if (polyscope::hasVolumeMesh(DELETE_PREVIEW)) {
                polyscope::removeStructure(polyscope::getVolumeMesh(DELETE_PREVIEW));
            }
            return;
        }

        const auto hexes = m_blocking->mesh_hexes(1);
        const int position = display_position(3, *m_hover_block);
        if (position < 0) {
            m_hover_block.reset();
            if (polyscope::hasVolumeMesh(DELETE_PREVIEW)) {
                polyscope::removeStructure(polyscope::getVolumeMesh(DELETE_PREVIEW));
            }
            return;
        }
        const auto index = static_cast<std::size_t>(position);
        if (index >= hexes.size()) {
            m_hover_block.reset();
            return;
        }
        // Its own 8 points rather than the whole blocking's: a structure carrying every corner would
        // draw a marker at each of them.
        const auto points = m_blocking->mesh_vertices(1);
        std::array<std::array<double, 3>, 8> corners{};
        std::array<double, 3> centre{0.0, 0.0, 0.0};
        for (std::size_t c = 0; c < 8; ++c) {
            corners[c] = points[static_cast<std::size_t>(hexes[index][c])];
            for (std::size_t k = 0; k < 3; ++k) {
                centre[k] += corners[c][k] / 8.0;
            }
        }

        // Grown very slightly about its own centre. At the block's own size the highlight sits
        // exactly on the faces it is meant to mark and z-fights with them — half its surface
        // winning, half losing, which reads on screen as no highlight rather than as one. A shell
        // just outside is unambiguous, and small enough not to look like a different block.
        const double scale = m_config.delete_highlight_scale;
        std::vector<std::array<double, 3>> shell;
        shell.reserve(8);
        for (const auto &corner : corners) {
            shell.push_back({centre[0] + (corner[0] - centre[0]) * scale,
                             centre[1] + (corner[1] - centre[1]) * scale,
                             centre[2] + (corner[2] - centre[2]) * scale});
        }

        auto *preview = polyscope::registerHexMesh(
            DELETE_PREVIEW, shell, std::vector<std::array<int, 8>>{{0, 1, 2, 3, 4, 5, 6, 7}});
        preview->setColor(to_glm(m_config.delete_highlight_color));
        preview->setEdgeWidth(1.0);
    }

    void BiyApp::handle_delete() {
        if (!m_blocking || m_mode != MouseMode::Delete) return;
        ImGuiIO &io = ImGui::GetIO();

        const bool by_key = !io.WantCaptureKeyboard && ImGui::IsKeyPressed(ImGuiKey_Space);
        const bool by_click = !io.WantCaptureMouse && io.MouseClicked[0];
        if (by_click || by_key) {
            // Acts on the block the highlight is showing, without re-testing what is under the
            // cursor — see handle_cut() for why that ordering matters on a trackpad.
            perform_delete(by_key ? "space" : "click");
            return;
        }
        if (io.WantCaptureMouse) return;

        const glm::vec2 mouse{io.MousePos.x, io.MousePos.y};
        if (mouse != m_last_delete_mouse) {
            m_last_delete_mouse = mouse;
            const auto previous = m_hover_block;
            update_delete_hover(mouse);
            refresh_delete_preview();
            // Said out loud as the cursor moves, not only once a block is gone: a highlight on its
            // own leaves "am I even on one?" unanswered precisely when the answer is no.
            if (m_hover_block != previous) {
                m_status = m_hover_block ? "Block " + std::to_string(*m_hover_block) + " — click or Space to delete"
                                         : "No block under the cursor";
            }
        }
    }

    void BiyApp::perform_delete(const char *trigger) {
        const auto report = [trigger](const std::string &message) {
            std::cout << "biy [delete by " << trigger << "]: " << message << std::endl;
        };

        if (!m_hover_block) {
            m_status = "Nothing under the cursor to delete";
            report("no block under the cursor. Point at one — it lights up when you are on it — then "
                   "click or press space. If \"blocks\" is unticked in the Scene panel there is "
                   "nothing to point at.");
            return;
        }

        const int block_id = *m_hover_block;
        const std::size_t before = m_blocking->nb_cells(3);
        m_blocking->delete_block(block_id);
        const std::size_t after = m_blocking->nb_cells(3);

        // The highlight describes a block that no longer exists, so it goes before anything is drawn
        // again.
        m_hover_block.reset();
        m_last_delete_mouse = glm::vec2(-1.0f, -1.0f);
        refresh_delete_preview();
        refresh_view();

        m_status = "Deleted block " + std::to_string(block_id) + " — blocks " + std::to_string(before) + " to " +
                   std::to_string(after);
        report("deleted block " + std::to_string(block_id) + ", blocks " + std::to_string(before) + " -> " +
               std::to_string(after) + "." + bend_report(*m_blocking));
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
