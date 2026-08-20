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
        /** @brief Color of the corner currently being dragged. White by default, so it can't be
         * mistaken for any of the classification colors below. */
        std::array<float, 3> corner_highlight_color{1.0f, 1.0f, 1.0f};

        /** @brief Color of a corner that isn't classified onto the model yet. */
        std::array<float, 3> corner_color_unclassified{0.6f, 0.2f, 0.85f};
        /** @brief Color of a corner classified on a model vertex. */
        std::array<float, 3> corner_color_on_vertex{1.0f, 0.9f, 0.1f};
        /** @brief Color of a corner classified on a model curve. */
        std::array<float, 3> corner_color_on_curve{0.9f, 0.15f, 0.15f};
        /** @brief Color of a corner classified on a model surface. */
        std::array<float, 3> corner_color_on_surface{0.15f, 0.4f, 0.95f};
        /** @brief Color of a corner classified on a model volume. Not one of the four states biy
         * was asked for, but classify() can produce it, so it gets its own color rather than being
         * silently drawn as unclassified. */
        std::array<float, 3> corner_color_on_volume{0.2f, 0.75f, 0.3f};

        /** @brief Whether the block structure's own edges are drawn at startup. */
        bool show_block_edges = true;
        /** @brief Thickness of those edges. */
        double block_edge_radius = 0.003;
        /** @brief Color of those edges. */
        std::array<float, 3> block_edge_color{0.15f, 0.15f, 0.15f};

        /**
         * @brief Picks the color a corner should be drawn in.
         * @param classification_dim Dimension the corner is classified on, or -1 if unclassified —
         *        as returned by `BlockingFacade::node_classification_dims()`.
         * @return The matching color.
         */
        [[nodiscard]] const std::array<float, 3> &color_for(int classification_dim) const;

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
