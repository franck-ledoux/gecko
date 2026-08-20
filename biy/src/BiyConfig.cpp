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
            read_color(j, "corner_color", config.corner_color);
            read_color(j, "corner_highlight_color", config.corner_highlight_color);
            message = "Loaded " + path;
        } catch (const nlohmann::json::exception &e) {
            // A broken config shouldn't stop biy from starting: say so and carry on with defaults.
            message = "Ignoring malformed " + path + " (" + e.what() + ")";
            return BiyConfig{};
        }

        return config;
    }

} // namespace gecko::biy
