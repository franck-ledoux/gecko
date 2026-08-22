# Python Bindings

Gecko ships a `gecko` Python extension module built with [pybind11](https://pybind11.readthedocs.io/).
Rather than binding every internal C++ class, it exposes a small number of **façade manager
classes** — `GeomModel` and `Blocking` today — in the spirit of
[gmsh's Python API](https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-API): every façade method takes
and returns only primitive types (`str`, `int`, `float`, and lists/tuples of those). Internal C++
types (`Point3d`, CGAL dart/attribute handles, `UnstructuredMesh`, ...) never cross into Python;
faces and blocks are identified to Python solely by plain `int` ids assigned by the façade.

!!! note
    This page, and the examples on it, are exercised by [`python/tst/`](https://github.com/franck-ledoux/gecko/tree/main/python/tst),
    which CI runs on every push. If this page and the actual behavior ever drift apart, those tests
    (and therefore CI) fail — so what's documented here is always what the code actually does.

## Building

The Python module is opt-in (it needs a Python development environment — headers included, not
just an interpreter):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DGECKO_BUILD_PYTHON=ON
cmake --build build --target gecko
```

This produces a `gecko` extension module (e.g. `gecko.cpython-312-x86_64-linux-gnu.so`) under
`build/python/`.

## Using it

Point `PYTHONPATH` at the directory containing the built module, then import it normally:

```bash
PYTHONPATH=build/python python3 -c "import gecko; print(gecko.hello())"
```

```
Hello from gecko!
```

### GeomModel

`GeomModel` loads a geometric model from a Gmsh MSH file (triangles for a boundary representation,
tetrahedra to also reconstruct volumes) and lets you inspect its entities and physical groups:

```python
import gecko

model = gecko.GeomModel("my_model.msh")
print(model.nb_vertices(), model.nb_curves(), model.nb_surfaces(), model.nb_volumes())
print(model.surface_tags())          # e.g. [1, 2, 3]

for group_id in model.group_ids():
    print(model.group_name(group_id), model.group_dim(group_id))
    print(model.group_entities(group_id))  # [(dim, entity_tag), ...]
```

### Blocking

`Blocking` builds a structured quad/hex blocking against a `GeomModel`. Block corners are passed
as lists of `(x, y, z)` tuples; created faces/blocks are identified by the `int` id `Blocking`
hands back. `degree=1` (the default) gives straight edges; `degree=3` gives cubic-Bezier edges that
`classify()` can bend onto curved geometry. Degrees 1 through 4 are available.

```python
import gecko

model = gecko.GeomModel("my_model.msh")
blocking = gecko.Blocking(model)  # model must stay alive for as long as blocking is used

face_a = blocking.create_quad_block([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)])
face_b = blocking.create_quad_block([(1, 0, 0), (2, 0, 0), (2, 1, 0), (1, 1, 0)])
blocking.build_connectivity()  # sews the shared edge between face_a and face_b

blocking.classify(tol_vertex=1e-6, tol_curve=1e-3, tol_surface=1e-2)  # snap onto the model

# What each cell ended up on: -1 free, 0 vertex, 1 curve, 2 surface, 3 volume. Edges and faces are
# inferred from their own boundary, not from proximity — see the biy guide.
print(blocking.node_classification_dims(), blocking.edge_classification_dims())

# Snap one corner back onto the model after moving it, reclassifying only what touches it.
blocking.move_node(0, 0.01, 0.01, 0.0)
blocking.snap_node(0, tol_vertex=0.1)

print(blocking.nb_cells(2), blocking.is_valid_topology())

if blocking.can_delete_face(face_a):
    blocking.delete_face(face_a)

# Cutting: pick an edge (by its position in edge_vertices()/edge_segments() order) and a parameter
# along it. The cut runs through every edge parallel to that one, across every block they touch —
# sheet_edges() reports which, and sheet_cut_points() where, before anything is modified.
sheet = blocking.sheet_edges(0)                 # [] if the sheet cannot be cut evenly
preview = blocking.sheet_cut_points(0, 0.5)     # one point per edge of `sheet`, in the same order
blocking.cut_sheet(0, 0.5)                      # False, changing nothing, if the cut is impossible

# What each block encloses, from its own stored geometry — exact at 1 subdivision when its faces are
# planar, converging as it grows otherwise. Negative means that block came out inverted.
print(blocking.block_volumes(subdivisions=1))

# The halves keep the geometry they were cut out of exactly (De Casteljau subdivision, not a refit),
# so meshing after a cut traces the same points as before. Note classify() rebuilds faces and blocks
# from their boundaries and so discards that — cut first, classify after.

# The exported mesh's nodes also carry 2 POINT_DATA scalars, "classification_dim"/"classification_tag"
# (-1/-1 when unclassified) — the same dim/tag pair node_classification_dims() reports, for every
# node the subdivision generates, not just the original block corners.
blocking.write_vtk(subdivisions=4, path="blocking.vtk")
```

## Running the tests

The Python test suite (`python/tst/`, using [pytest](https://pytest.org)) is registered as a
regular CTest test (`test_python`) alongside every C++ suite, once `GECKO_BUILD_TESTS` (the
default) and `GECKO_BUILD_PYTHON` are both on:

```bash
pip install -r python/requirements-test.txt
cmake --build build
cd build && ctest -R test_python --output-on-failure
```

## Code coverage

Because `gecko`'s C++ source (`python/src/bindings.cpp`) is compiled by the very same CMake build
as the rest of the library, it picks up the project's `--coverage` instrumentation automatically
whenever `GECKO_CODE_COVERAGE`/`GECKO_CODE_COVERAGE_HTML_REPORT` is on (see
[Building & Testing](../developer-guide/testing.md#code-coverage-html-report)) — and because
`test_python` is a regular CTest test, running it exercises that instrumented code exactly like any
C++ test does. In other words: **Python tests contribute to the same C++ coverage report as the
C++ tests**, covering both the binding glue code and whatever C++ library code it calls into — no
separate coverage tooling needed for the Python side.
