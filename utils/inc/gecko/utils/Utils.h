#pragma once
#include <cassert>

namespace gecko {
    /**
     * @brief Signals that a code path is logically unreachable and must never execute.
     *
     * Asserts in Debug builds, then hints the compiler that the path is dead so it
     * can optimize accordingly (e.g. after an exhaustive switch on an enum).
     */
    [[noreturn]] inline void unreachable() noexcept {
        assert(false && "Unreachable code path executed!");
#if defined(__GNUC__) || defined(__clang__)
        __builtin_unreachable();
#elif defined(_MSC_VER)
        __assume(0);
#endif
    }
} // namespace gecko
