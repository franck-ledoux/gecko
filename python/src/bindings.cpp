// gecko Python bindings — currently a minimal "hello world" scaffold: it exists to prove the
// build/test/coverage/doc infrastructure works end-to-end, not to expose the C++ API yet. See
// docs/user-guide/python.md.
#include <string>

#include <pybind11/pybind11.h>

namespace {
    std::string hello() { return "Hello from gecko!"; }
} // namespace

PYBIND11_MODULE(gecko, m) {
    m.doc() = "Gecko Python bindings (early scaffolding, see issue #25).";
    m.def("hello", &hello, "Returns a greeting string, confirming the gecko extension loads and runs.");
}
