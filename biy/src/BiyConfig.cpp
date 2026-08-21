#include "BiyConfig.h"

#include <filesystem>
#include <fstream>

#include <nlohmann/json.hpp>

namespace gecko::biy {

    namespace {
        /** @brief Reads one numeric key into @p ATarget, leaving it untouched if absent. */
        void read_double(const nlohmann::json &j, const char *key, double &target) {
            if (j.contains(key) && j.at(key).is_number()) target = j.at(key).get<double>();
        }

        /** @brief Reads one boolean key into @p ATarget, leaving it untouched if absent. */
        void read_bool(const nlohmann::json &j, const char *key, bool &target) {
            if (j.contains(key) && j.at(key).is_boolean()) target = j.at(key).get<bool>();
        }

        /** @brief Reads one `[r,g,b]` key into @p ATarget, leaving it untouched if absent. */
        void read_color(const nlohmann::json &j, const char *key, std::array<float, 3> &target) {
            if (!j.contains(key)) return;
            const auto &value = j.at(key);
            if (!value.is_array() || value.size() != 3) return;
            for (std::size_t i = 0; i < 3; ++i) {
                if (value[i].is_number()) target[i] = value[i].get<float>();
            }
        }
    } // namespace

    const std::array<float, 3> &BiyConfig::color_for(int classification_dim) const {
        switch (classification_dim) {
            case 0:
                return corner_color_on_vertex;
            case 1:
                return corner_color_on_curve;
            case 2:
                return corner_color_on_surface;
            case 3:
                return corner_color_on_volume;
            default:
                return corner_color_unclassified;
        }
    }

    BiyConfig BiyConfig::load(const std::string &path, std::string &message) {
        BiyConfig config;

        if (!std::filesystem::exists(path)) {
            message = "No " + path + " found, using defaults";
            return config;
        }

        std::ifstream file(path);
        if (!file) {
            message = "Could not open " + path + ", using defaults";
            return config;
        }

        try {
            const nlohmann::json j = nlohmann::json::parse(file);
            read_double(j, "corner_radius", config.corner_radius);
            read_double(j, "corner_highlight_radius", config.corner_highlight_radius);
            read_color(j, "corner_highlight_color", config.corner_highlight_color);
            read_color(j, "corner_color_unclassified", config.corner_color_unclassified);
            read_color(j, "corner_color_on_vertex", config.corner_color_on_vertex);
            read_color(j, "corner_color_on_curve", config.corner_color_on_curve);
            read_color(j, "corner_color_on_surface", config.corner_color_on_surface);
            read_color(j, "corner_color_on_volume", config.corner_color_on_volume);
            read_double(j, "model_transparency", config.model_transparency);
            read_double(j, "panel_width", config.panel_width);
            read_double(j, "gizmo_size", config.gizmo_size);
            read_double(j, "gizmo_dot_radius", config.gizmo_dot_radius);
            read_color(j, "gizmo_color_x", config.gizmo_color_x);
            read_color(j, "gizmo_color_y", config.gizmo_color_y);
            read_color(j, "gizmo_color_z", config.gizmo_color_z);
            read_color(j, "geometry_volume_color", config.geometry_volume_color);
            read_color(j, "geometry_surface_color", config.geometry_surface_color);
            read_color(j, "geometry_curve_color", config.geometry_curve_color);
            read_double(j, "geometry_curve_radius", config.geometry_curve_radius);
            read_color(j, "geometry_point_color", config.geometry_point_color);
            read_double(j, "geometry_point_radius", config.geometry_point_radius);
            read_color(j, "mesh_color", config.mesh_color);
            read_bool(j, "show_geometry_volumes", config.show_geometry_volumes);
            read_bool(j, "show_geometry_surfaces", config.show_geometry_surfaces);
            read_bool(j, "show_geometry_curves", config.show_geometry_curves);
            read_bool(j, "show_geometry_points", config.show_geometry_points);
            read_bool(j, "show_mesh", config.show_mesh);
            read_bool(j, "show_edge_control", config.show_edge_control);
            read_bool(j, "show_face_control", config.show_face_control);
            read_bool(j, "show_block_control", config.show_block_control);
            read_double(j, "control_point_radius", config.control_point_radius);
            read_double(j, "control_polygon_radius", config.control_polygon_radius);
            read_color(j, "edge_control_color", config.edge_control_color);
            read_color(j, "face_control_color", config.face_control_color);
            read_color(j, "block_control_color", config.block_control_color);
            read_double(j, "tol_vertex", config.tol_vertex);
            read_double(j, "tol_curve", config.tol_curve);
            read_double(j, "tol_surface", config.tol_surface);
            read_bool(j, "show_block_edges", config.show_block_edges);
            read_double(j, "block_edge_radius", config.block_edge_radius);
            read_color(j, "block_edge_color", config.block_edge_color);
            message = "Loaded " + path;
        } catch (const nlohmann::json::exception &e) {
            // A broken config shouldn't stop biy from starting: say so and carry on with defaults.
            message = "Ignoring malformed " + path + " (" + e.what() + ")";
            return BiyConfig{};
        }

        return config;
    }

} // namespace gecko::biy
