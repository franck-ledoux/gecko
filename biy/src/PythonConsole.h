#pragma once

#include <atomic>
#include <thread>

namespace gecko::biy {

    class BiyApp;

    /**
     * @class PythonConsole
     * @brief A Python REPL over stdin, running on its own thread, with biy's live `model` and
     * `blocking` objects already bound in its namespace.
     *
     * The interpreter is embedded (not a subprocess): statements typed here mutate the very same
     * C++ facade objects the 3D view displays, with no serialization in between. Each statement
     * runs while holding the app's mutex, so it never races the render thread; the view is
     * refreshed afterwards so console-driven changes appear immediately.
     */
    class PythonConsole {
    public:
        /**
         * @brief Constructor. Starts the REPL thread.
         * @param app The application whose state to expose; must outlive this console.
         */
        explicit PythonConsole(BiyApp &app);

        /** @brief Destructor. Signals the REPL thread to stop and joins it. */
        ~PythonConsole();

        PythonConsole(const PythonConsole &) = delete;
        PythonConsole &operator=(const PythonConsole &) = delete;

        /** @brief Asks the REPL thread to stop at its next opportunity. */
        void stop() { m_running = false; }

        /** @brief Checks whether the user ended the session from the console (Ctrl-D / `exit()`).
         * @return true while the REPL is still running. */
        [[nodiscard]] bool running() const { return m_running; }

    private:
        void repl(BiyApp &app);

        std::atomic<bool> m_running{true};
        std::thread m_thread;
    };

} // namespace gecko::biy
