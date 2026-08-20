#include "PythonConsole.h"

#include <iostream>
#include <mutex>
#include <string>

#include <pybind11/embed.h>

#include "BiyApp.h"
#include "PyGeckoBindings.h"

namespace py = pybind11;

// The same bindings the standalone `gecko` extension module exposes, registered here as a built-in
// of the embedded interpreter — so `import gecko` works inside biy without locating a built .so,
// and yields the exact same GeomModel/Blocking classes.
PYBIND11_EMBEDDED_MODULE(gecko, m) { gecko::python::register_gecko_module(m); }

namespace gecko::biy {

    PythonConsole::PythonConsole(BiyApp &app) : m_thread([this, &app] { repl(app); }) {}

    PythonConsole::~PythonConsole() {
        m_running = false;
        if (m_thread.joinable()) m_thread.join();
    }

    void PythonConsole::repl(BiyApp &app) {
        // The interpreter lives entirely on this thread: the render thread never touches Python, so
        // no GIL juggling is needed beyond what pybind11 does for us here.
        const py::scoped_interpreter interpreter;

        py::dict globals = py::module_::import("__main__").attr("__dict__");
        {
            const std::lock_guard<std::mutex> lock(app.mutex());
            globals["gecko"] = py::module_::import("gecko");
            // Bound by reference: these are the very objects the 3D view is displaying, not copies.
            globals["model"] = py::cast(&app.model(), py::return_value_policy::reference);
            globals["blocking"] = py::cast(&app.blocking(), py::return_value_policy::reference);
        }

        std::cout << "biy Python console — `model` and `blocking` are live; Ctrl-D to leave.\n"
                  << ">>> " << std::flush;

        std::string line;
        while (m_running && std::getline(std::cin, line)) {
            if (!m_running) break;
            if (!line.empty()) {
                const std::lock_guard<std::mutex> lock(app.mutex());
                try {
                    // Evaluate as an expression first so the REPL echoes results like Python's own
                    // does; statements (assignments, loops) aren't expressions, so fall back to exec.
                    try {
                        const py::object result = py::eval(line, globals);
                        if (!result.is_none()) std::cout << py::str(result).cast<std::string>() << "\n";
                    } catch (const py::error_already_set &) {
                        py::exec(line, globals);
                    }
                    app.refresh_view();
                } catch (const py::error_already_set &e) {
                    std::cout << e.what() << "\n";
                } catch (const std::exception &e) {
                    std::cout << "error: " << e.what() << "\n";
                }
            }
            std::cout << ">>> " << std::flush;
        }

        m_running = false;
        std::cout << "\nConsole closed.\n" << std::flush;
    }

} // namespace gecko::biy
