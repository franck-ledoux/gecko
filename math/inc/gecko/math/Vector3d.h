#pragma once

#include <cmath>
#include <cstddef>
#include <utility>
#include <gecko/utils/Types.h>
#include <gecko/math/Point3d.h>

namespace gecko {
    // =========================================================================
    // CLASS Vector3d
    // =========================================================================
    /**
     * @class Vector3d
     * @brief Represents a vector in 3D space, defined by its (x, y, z) components.
     */
    class Vector3d {
    public:
        /** @brief Default constructor. Builds the null vector (0, 0, 0). */
        constexpr Vector3d() noexcept = default;

        /**
         * @brief Constructor from explicit components.
         * @param x_val X component.
         * @param y_val Y component.
         * @param z_val Z component (default: 0.0, useful for 2D vectors).
         */
        constexpr Vector3d(Float x_val, Float y_val, Float z_val = 0.0) noexcept : m_x(x_val), m_y(y_val), m_z(z_val) {}

        /**
         * @brief Constructs the vector going from a point p1 to a point p2 (p2 - p1).
         * @param p1 Start point.
         * @param p2 End point.
         */
        constexpr Vector3d(const Point3d &p1, const Point3d &p2) noexcept
            : m_x(p2.x() - p1.x()), m_y(p2.y() - p1.y()), m_z(p2.z() - p1.z()) {}

        /** @brief Gets the X component. @return The X component. */
        [[nodiscard]] constexpr Float x() const noexcept { return m_x; }
        /** @brief Gets the Y component. @return The Y component. */
        [[nodiscard]] constexpr Float y() const noexcept { return m_y; }
        /** @brief Gets the Z component. @return The Z component. */
        [[nodiscard]] constexpr Float z() const noexcept { return m_z; }

        /** @brief Gets a mutable reference to the X component. @return Reference to the X component. */
        [[nodiscard]] constexpr Float &x() noexcept { return m_x; }
        /** @brief Gets a mutable reference to the Y component. @return Reference to the Y component. */
        [[nodiscard]] constexpr Float &y() noexcept { return m_y; }
        /** @brief Gets a mutable reference to the Z component. @return Reference to the Z component. */
        [[nodiscard]] constexpr Float &z() noexcept { return m_z; }

        /**
         * @brief Accesses a component by index (0 = x, 1 = y, 2 = z).
         * @param i Index of the requested component.
         * @return The value of the requested component.
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
                    std::unreachable();
            }
        }

        /**
         * @brief Accesses a component by index (0 = x, 1 = y, 2 = z).
         * @param i Index of the requested component.
         * @return Reference to the requested component.
         * @pre @p i must be in [0, 2].
         */
        [[nodiscard]] constexpr Float &operator[](std::size_t i) noexcept {
            switch (i) {
                case 0:
                    return m_x;
                case 1:
                    return m_y;
                case 2:
                    return m_z;
                default:
                    std::unreachable();
            }
        }

        /**
         * @brief Equality comparison (C++20 defaulted, compares all components).
         * @return true if both vectors have the same components, false otherwise.
         */
        bool operator==(const Vector3d &) const = default;

        /**
         * @brief Unary negation.
         * @return The opposite vector (-x, -y, -z).
         */
        constexpr Vector3d operator-() const noexcept { return {-m_x, -m_y, -m_z}; }

        /**
         * @brief Adds another vector to this one in place.
         * @param v Vector to add.
         * @return Reference to this vector after the addition.
         */
        constexpr Vector3d &operator+=(const Vector3d &v) noexcept {
            m_x += v.m_x;
            m_y += v.m_y;
            m_z += v.m_z;
            return *this;
        }

        /**
         * @brief Subtracts another vector from this one in place.
         * @param v Vector to subtract.
         * @return Reference to this vector after the subtraction.
         */
        constexpr Vector3d &operator-=(const Vector3d &v) noexcept {
            m_x -= v.m_x;
            m_y -= v.m_y;
            m_z -= v.m_z;
            return *this;
        }

        /**
         * @brief Scales this vector in place by a scalar factor.
         * @param scalar Scaling factor.
         * @return Reference to this vector after scaling.
         */
        constexpr Vector3d &operator*=(Float scalar) noexcept {
            m_x *= scalar;
            m_y *= scalar;
            m_z *= scalar;
            return *this;
        }

        /**
         * @brief Divides this vector in place by a scalar factor.
         * @param scalar Divisor.
         * @return Reference to this vector after the division.
         * @pre @p scalar must not be 0.
         */
        constexpr Vector3d &operator/=(Float scalar) noexcept {
            m_x /= scalar;
            m_y /= scalar;
            m_z /= scalar;
            return *this;
        }

        /**
         * @brief Computes the dot product between this vector and another one.
         * @param v The other vector.
         * @return The scalar dot product.
         */
        [[nodiscard]] constexpr Float dot(const Vector3d &v) const noexcept {
            return m_x * v.m_x + m_y * v.m_y + m_z * v.m_z;
        }

        /**
         * @brief Computes the cross product between this vector and another one.
         * @param v The other vector.
         * @return The resulting vector, orthogonal to both this vector and @p v.
         */
        [[nodiscard]] constexpr Vector3d cross(const Vector3d &v) const noexcept {
            return {m_y * v.m_z - m_z * v.m_y, m_z * v.m_x - m_x * v.m_z, m_x * v.m_y - m_y * v.m_x};
        }

        /**
         * @brief Computes the squared Euclidean norm of this vector.
         * @return The squared norm (avoids the cost of a square root).
         */
        [[nodiscard]] constexpr Float norm_sq() const noexcept { return dot(*this); }

        /**
         * @brief Computes the Euclidean norm (length) of this vector.
         * @return The norm of the vector.
         */
        [[nodiscard]] Float norm() const noexcept { return std::sqrt(norm_sq()); }

        /**
         * @brief Computes a normalized (unit length) copy of this vector.
         * @return The normalized vector, or the null vector if this vector's norm is 0.
         */
        [[nodiscard]] Vector3d normalized() const noexcept {
            double n = norm();
            if (n > 0.0) {
                Vector3d res = *this;
                res /= n; // Utilise operator/= membre
                return res;
            }
            return Vector3d{};
        }

    private:
        Float m_x{0.0}, m_y{0.0}, m_z{0.0};
    };

    // =========================================================================

    // See Point3d::operator+=() in Point3d.h for the documentation of this out-of-line definition.
    constexpr Point3d &Point3d::operator+=(const Vector3d &v) noexcept {
        m_x += v.x();
        m_y += v.y();
        m_z += v.z();
        return *this;
    }

    // See Point3d::operator-=() in Point3d.h for the documentation of this out-of-line definition.
    constexpr Point3d &Point3d::operator-=(const Vector3d &v) noexcept {
        m_x -= v.x();
        m_y -= v.y();
        m_z -= v.z();
        return *this;
    }

    /**
     * @brief Computes the vector going from p2 to p1 (Point3d - Point3d = Vector3d).
     * @param p1 First point.
     * @param p2 Second point.
     * @return The vector p1 - p2.
     */
    [[nodiscard]] constexpr Vector3d operator-(const Point3d &p1, const Point3d &p2) noexcept {
        return {p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z()};
    }

    /**
     * @brief Translates a point by a vector (Point3d + Vector3d = Point3d).
     * @param p Point to translate.
     * @param v Translation vector.
     * @return The translated point.
     */
    [[nodiscard]] constexpr Point3d operator+(Point3d p, const Vector3d &v) noexcept {
        p += v;
        return p;
    }

    /**
     * @brief Translates a point by a vector (Vector3d + Point3d = Point3d).
     * @param v Translation vector.
     * @param p Point to translate.
     * @return The translated point.
     */
    [[nodiscard]] constexpr Point3d operator+(const Vector3d &v, Point3d p) noexcept {
        p += v;
        return p;
    }

    /**
     * @brief Translates a point by the opposite of a vector (Point3d - Vector3d = Point3d).
     * @param p Point to translate.
     * @param v Translation vector to subtract.
     * @return The translated point.
     */
    [[nodiscard]] constexpr Point3d operator-(Point3d p, const Vector3d &v) noexcept {
        p -= v;
        return p;
    }

    /**
     * @brief Adds two vectors (Vector3d + Vector3d = Vector3d).
     * @param v1 First vector.
     * @param v2 Second vector.
     * @return The sum vector v1 + v2.
     */
    [[nodiscard]] constexpr Vector3d operator+(Vector3d v1, const Vector3d &v2) noexcept {
        v1 += v2;
        return v1;
    }

    /**
     * @brief Subtracts two vectors (Vector3d - Vector3d = Vector3d).
     * @param v1 First vector.
     * @param v2 Second vector.
     * @return The difference vector v1 - v2.
     */
    [[nodiscard]] constexpr Vector3d operator-(Vector3d v1, const Vector3d &v2) noexcept {
        v1 -= v2;
        return v1;
    }

    /**
     * @brief Scales a vector by a scalar factor (Vector3d * Float).
     * @param v Vector to scale.
     * @param scalar Scaling factor.
     * @return The scaled vector.
     */
    [[nodiscard]] constexpr Vector3d operator*(Vector3d v, Float scalar) noexcept {
        v *= scalar;
        return v;
    }

    /**
     * @brief Divides a vector by a scalar factor (Vector3d / Float).
     * @param v Vector to divide.
     * @param scalar Divisor.
     * @return The scaled vector.
     * @pre @p scalar must not be 0.
     */
    [[nodiscard]] constexpr Vector3d operator/(Vector3d v, Float scalar) noexcept {
        v /= scalar;
        return v;
    }

    /**
     * @brief Scales a vector by a scalar factor (Float * Vector3d).
     * @param scalar Scaling factor.
     * @param v Vector to scale.
     * @return The scaled vector.
     */
    [[nodiscard]] constexpr Vector3d operator*(Float scalar, Vector3d v) noexcept {
        v *= scalar;
        return v;
    }

    /**
     * @brief Computes the dot product between two vectors.
     * @param u First vector.
     * @param v Second vector.
     * @return The scalar dot product.
     * @see Vector3d::dot()
     */
    [[nodiscard]] constexpr Float dot(const Vector3d &u, const Vector3d &v) noexcept { return u.dot(v); }

    /**
     * @brief Computes the cross product between two vectors.
     * @param u First vector.
     * @param v Second vector.
     * @return The resulting vector, orthogonal to both @p u and @p v.
     * @see Vector3d::cross()
     */
    [[nodiscard]] constexpr Vector3d cross(const Vector3d &u, const Vector3d &v) noexcept { return u.cross(v); }
} // namespace gecko
