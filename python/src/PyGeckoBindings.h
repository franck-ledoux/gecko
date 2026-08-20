#pragma once

#include <pybind11/pybind11.h>

namespace gecko::python {

    /**
     * @brief Registers gecko's façade classes (GeomModel, Blocking) into a Python module.
     *
     * Shared by the standalone `gecko` extension module (python/src/bindings.cpp) and by biy's
     * embedded interpreter (biy/src/PythonConsole.cpp), so both expose the exact same API from a
     * single definition.
     * @param m The module to register into.
     */
    void register_gecko_module(pybind11::module_ m);

} // namespace gecko::python
