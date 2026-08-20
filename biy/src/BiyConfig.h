#pragma once

#include <array>
#include <string>

namespace gecko::biy {

    /**
     * @struct BiyConfig
     * @brief biy's tunable display settings, read from a `biy_config.json` file at startup.
     *
     * Every field has a usable default, so the file is entirely optional and may set only the keys
     * it cares about. Radii are Polyscope "relative" values (a fraction of the scene's bounding
     * box), so they stay sensible whatever the model's units are; colors are RGB in [0,1].
     */
    struct BiyConfig {
        /** @brief Radius of a block corner at rest. */
        double corner_radius = 0.01;
        /** @brief Radius of the corner currently being dragged. */
        double corner_highlight_radius = 0.02;
        /** @brief Color of a block corner at rest. */
        std::array<float, 3> corner_color{0.85f, 0.15f, 0.75f};
        /** @brief Color of the corner currently being dragged. */
        std::array<float, 3> corner_highlight_color{1.0f, 0.85f, 0.1f};

        /**
         * @brief Loads settings from a JSON file, falling back to the defaults above for the file's
         * missing (or malformed) keys.
         * @param path Path to the JSON file; a missing file is not an error.
         * @param message Set to a human-readable description of what was loaded, or of why the file
         *        was ignored — shown in biy's panel.
         * @return The resulting configuration.
         */
        static BiyConfig load(const std::string &path, std::string &message);
    };

} // namespace gecko::biy
