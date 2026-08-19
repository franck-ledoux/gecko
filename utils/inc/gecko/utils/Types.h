#pragma once
#include <compare>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace gecko {
    using DefaultIndexType = std::uint32_t; // Or std::size_t
    /**
     * @class StrongId
     * @brief Template structure to get a strong typing for cells and avoid to mix node, edge and face tags.
     * @tparam Tag Empty tag type used to distinguish StrongId instantiations at compile-time.
     * @tparam T Underlying integral type used to store the id value.
     */
    template<typename Tag, typename T = DefaultIndexType>
    struct StrongId {
        /** @brief Underlying integral type used to store the id value. */
        using value_type = T;
        /** @brief Sentinel value used to represent an invalid id. */
        static constexpr T invalid_value = std::numeric_limits<T>::max();

        /** @brief Underlying id value, defaults to #invalid_value. */
        T value{invalid_value};

        /**
         * @brief Sentinel invalid StrongId instance.
         * @note Declared here without an initializer and defined just after the class: a static
         *       data member of the class's own type cannot be initialized in-class, because the
         *       class is still an incomplete type at that point (and `constexpr` requires an
         *       initializer right where it's declared).
         */
        static const StrongId invalid_id;

        /** @brief Default constructor, builds an invalid id. */
        constexpr StrongId() noexcept : value(invalid_value) {}
        /**
         * @brief Constructs an id from its underlying value.
         * @param val The underlying value to wrap.
         */
        constexpr explicit StrongId(T val) noexcept : value(val) {}

        /** @brief Resets the id to the invalid sentinel value. */
        constexpr void reset() noexcept { value = invalid_value; }

        /**
         * @brief Checks whether the id holds a value different from the invalid sentinel.
         * @return true if the id is valid, false otherwise.
         */
        [[nodiscard]] constexpr bool is_valid() const noexcept { return value != invalid_value; }

        /**
         * @brief Explicit conversion to the underlying value type (e.g. to index a vector).
         * @return The underlying value.
         */
        constexpr explicit operator T() const noexcept { return value; }

        // ---------------------------------------------------------------------
        // 2. Comparaisons (C++20)
        // ---------------------------------------------------------------------
        /**
         * @brief Three-way comparison, automatically provides ==, !=, <, <=, >, >=
         * between two StrongId of the same type.
         * @return The ordering between this id and the other one.
         */
        constexpr auto operator<=>(const StrongId &) const = default;

        // ---------------------------------------------------------------------
        // 3. Incrémentation & Décrémentation
        // ---------------------------------------------------------------------
        /**
         * @brief Pre-increment operator (++id).
         * @return Reference to this id after incrementing its value.
         */
        constexpr StrongId &operator++() noexcept {
            ++value;
            return *this;
        }

        /**
         * @brief Post-increment operator (id++).
         * @return A copy of this id before incrementing its value.
         */
        constexpr StrongId operator++(int) noexcept {
            StrongId tmp = *this;
            ++value;
            return tmp;
        }

        /**
         * @brief Pre-decrement operator (--id).
         * @return Reference to this id after decrementing its value.
         */
        constexpr StrongId &operator--() noexcept {
            --value;
            return *this;
        }

        /**
         * @brief Post-decrement operator (id--).
         * @return A copy of this id before decrementing its value.
         */
        constexpr StrongId operator--(int) noexcept {
            StrongId tmp = *this;
            --value;
            return tmp;
        }

        // ---------------------------------------------------------------------
        // 4. Opérateurs d'affectation arithmétique (Optionnel mais très pratique)
        // ---------------------------------------------------------------------
        /**
         * @brief Adds an offset to the id's underlying value in place.
         * @tparam OffT Type of the offset.
         * @param offset Offset to add.
         * @return Reference to this id after the addition.
         */
        template<typename OffT>
        constexpr StrongId &operator+=(OffT offset) noexcept {
            value = static_cast<T>(static_cast<OffT>(value) + offset);
            return *this;
        }

        /**
         * @brief Subtracts an offset from the id's underlying value in place.
         * @tparam OffT Type of the offset.
         * @param offset Offset to subtract.
         * @return Reference to this id after the subtraction.
         */
        template<typename OffT>
        constexpr StrongId &operator-=(OffT offset) noexcept {
            value = static_cast<T>(static_cast<OffT>(value) - offset);
            return *this;
        }
    };

    /** @brief Out-of-line definition of StrongId::invalid_id, now that StrongId is complete. */
    template<typename Tag, typename T>
    inline constexpr StrongId<Tag, T> StrongId<Tag, T>::invalid_id{};

    /**
     * @brief Adds an offset to a StrongId (id + offset).
     * @tparam TagT Tag type of the StrongId.
     * @tparam T Underlying value type of the StrongId.
     * @tparam OffT Type of the offset.
     * @param id The id to offset.
     * @param offset Offset to add.
     * @return A new StrongId whose value is id + offset.
     */
    template<typename TagT, typename T, typename OffT>
    constexpr StrongId<TagT, T> operator+(StrongId<TagT, T> id, OffT offset) noexcept {
        id += offset;
        return id;
    }

    /**
     * @brief Adds an offset to a StrongId (offset + id).
     * @tparam TagT Tag type of the StrongId.
     * @tparam T Underlying value type of the StrongId.
     * @tparam OffT Type of the offset.
     * @param offset Offset to add.
     * @param id The id to offset.
     * @return A new StrongId whose value is id + offset.
     */
    template<typename TagT, typename T, typename OffT>
    constexpr StrongId<TagT, T> operator+(OffT offset, StrongId<TagT, T> id) noexcept {
        id += offset;
        return id;
    }

    /**
     * @brief Subtracts an offset from a StrongId (id - offset).
     * @tparam TagT Tag type of the StrongId.
     * @tparam T Underlying value type of the StrongId.
     * @tparam OffT Type of the offset.
     * @param id The id to offset.
     * @param offset Offset to subtract.
     * @return A new StrongId whose value is id - offset.
     */
    template<typename TagT, typename T, typename OffT>
    constexpr StrongId<TagT, T> operator-(StrongId<TagT, T> id, OffT offset) noexcept {
        id -= offset;
        return id;
    }

    // Mesh type
    //Declaration of one structure per each cell type
    struct NodeTag {};
    struct EdgeTag {};
    struct FaceTag {};
    struct CellTag {};
    // By construction, the following type are distinct
    using NodeId = StrongId<NodeTag>;
    using EdgeId = StrongId<EdgeTag>;
    using FaceId = StrongId<FaceTag>;
    using CellId = StrongId<CellTag>;

    using SizeT = std::size_t;
    using Int = std::int32_t;
    using UInt = std::uint32_t;
    using Float = double;
} // namespace gecko
