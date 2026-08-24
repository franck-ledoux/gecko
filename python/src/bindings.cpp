// gecko Python bindings: a façade over gecko's C++ API, exposing only manager classes (GeomModel,
// Blocking) whose methods take/return primitive types (str/int/float and collections thereof) —
// never a raw C++ type (Point3d, CGAL handle, UnstructuredMesh, ...). See docs/user-guide/python.md.
//
// The bindings themselves live in PyGeckoBindings.cpp, shared with biy's embedded interpreter.
#include <pybind11/pybind11.h>

#include "PyGeckoBindings.h"

PYBIND11_MODULE(gecko, m) { gecko::app::register_gecko_module(m); }
