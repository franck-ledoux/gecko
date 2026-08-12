#pragma once
#include <cstddef>

#include <gecko/utils/Types.h>
#include <gecko/utils/Utils.h>

namespace gecko {

    class Vector3d;

    /**
     * @class Point3d
     * @brief Represents a point in 3D space, defined by its (x, y, z) coordinates.
     */
    class Point3d {

    public:
        /** @brief Default constructor. Builds the point (0, 0, 0). */
        constexpr Point3d() noexcept = default;

        /**
         * @brief Constructor from explicit coordinates.
         * @param x_val X coordinate.
         * @param y_val Y coordinate.
         * @param z_val Z coordinate (default: 0.0, useful for 2D points).
         */
        constexpr Point3d(Float x_val, Float y_val, Float z_val = 0.0) noexcept : m_x(x_val), m_y(y_val), m_z(z_val) {}

        /** @brief Gets the X coordinate. @return The X coordinate. */
        [[nodiscard]] constexpr Float x() const noexcept { return m_x; }
        /** @brief Gets the Y coordinate. @return The Y coordinate. */
        [[nodiscard]] constexpr Float y() const noexcept { return m_y; }
        /** @brief Gets the Z coordinate. @return The Z coordinate. */
        [[nodiscard]] constexpr Float z() const noexcept { return m_z; }

        /** @brief Gets a mutable reference to the X coordinate. @return Reference to the X coordinate. */
        [[nodiscard]] constexpr Float &x() noexcept { return m_x; }
        /** @brief Gets a mutable reference to the Y coordinate. @return Reference to the Y coordinate. */
        [[nodiscard]] constexpr Float &y() noexcept { return m_y; }
        /** @brief Gets a mutable reference to the Z coordinate. @return Reference to the Z coordinate. */
        [[nodiscard]] constexpr Float &z() noexcept { return m_z; }

        /**
         * @brief Accesses a coordinate by index (0 = x, 1 = y, 2 = z).
         * @param i Index of the requested coordinate.
         * @return The value of the requested coordinate.
         * @pre @p i must be in [0, 2].
         */
        [[nodiscard]] constexpr Float operator[](std::size_t i) const noexcept {
            switch (i) {
                case 0:
                    return m_x;
                case 1:
                    return m_y;
                case 2:
                    return m_z;
                default:
                    unreachable();
            }
        }

        /**
         * @brief Accesses a coordinate by index (0 = x, 1 = y, 2 = z).
         * @param i Index of the requested coordinate.
         * @return Reference to the requested coordinate.
         * @pre @p i must be in [0, 2].
         */
        [[nodiscard]] constexpr double &operator[](std::size_t i) noexcept {
            switch (i) {
                case 0:
                    return m_x;
                case 1:
                    return m_y;
                case 2:
                    return m_z;
                default:
                    unreachable();
            }
        }

        /**
         * @brief Equality comparison (C++20 defaulted, compares all coordinates).
         * @return true if both points have the same coordinates, false otherwise.
         */
        bool operator==(const Point3d &) const = default;

        /**
         * @brief Translates this point in place by a vector.
         * @param v Translation vector.
         * @return Reference to this point after translation.
         */
        constexpr Point3d &operator+=(const Vector3d &v) noexcept;

        /**
         * @brief Translates this point in place by the opposite of a vector.
         * @param v Translation vector to subtract.
         * @return Reference to this point after translation.
         */
        constexpr Point3d &operator-=(const Vector3d &v) noexcept;

    private:
        Float m_x{0.0}, m_y{0.0}, m_z{0.0};
    };
} // namespace gecko
