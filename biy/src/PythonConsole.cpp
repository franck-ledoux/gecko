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

        // Python's own REPL uses codeop to decide whether what's been typed so far is a complete
        // statement — reuse it rather than guessing, so `for`/`if`/`def` blocks and bracket
        // continuations behave here exactly as they do in `python3`.
        py::object compile_command = py::module_::import("codeop").attr("compile_command");
        py::object builtins_eval = py::module_::import("builtins").attr("eval");

        std::cout << "biy Python console — `model` and `blocking` are live; Ctrl-D to leave.\n"
                  << ">>> " << std::flush;

        std::string line;
        std::string pending;
        while (m_running && std::getline(std::cin, line)) {
            if (!m_running) break;

            // Joined without a trailing newline, the way code.InteractiveConsole does it: a
            // trailing newline tells codeop the input is finished, so a block header would be
            // rejected for having no body instead of asking for more lines.
            const bool continuing = !pending.empty();
            if (continuing) pending += "\n";
            pending += line;

            if (line.empty() && pending.find_first_not_of(" \t\r\n") == std::string::npos) {
                pending.clear();
                std::cout << ">>> " << std::flush;
                continue;
            }

            bool incomplete = false;
            {
                const std::lock_guard<std::mutex> lock(app.mutex());
                try {
                    // "single" mode makes Python print expression results for us, exactly like its
                    // own interactive prompt; a None result means "keep reading".
                    const py::object code = compile_command(pending, "<biy>", "single");
                    if (code.is_none()) {
                        incomplete = true;
                    } else {
                        builtins_eval(code, globals);
                        app.refresh_view();
                    }
                } catch (const py::error_already_set &e) {
                    std::cout << e.what() << "\n";
                } catch (const std::exception &e) {
                    std::cout << "error: " << e.what() << "\n";
                }
            }

            if (incomplete) {
                std::cout << "... " << std::flush;
                continue;
            }
            pending.clear();
            std::cout << ">>> " << std::flush;
        }

        m_running = false;
        std::cout << "\nConsole closed.\n" << std::flush;
    }

} // namespace gecko::biy
