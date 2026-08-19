# Gecko

GECKO stands for hi**G**h ord**E**r blo**C**K t**O**olkit. It provides data structures and
algorithms to easily modify a block structure used for discretizing a geometrical domain with
a hexahedral mesh.

This site is the **developer documentation** for Gecko: it covers how to build the project, how
to write and run its unit tests, how to measure code coverage, how to write and generate this
very documentation, and how code/docs formatting is enforced.

- Start with the [Building & Testing](developer-guide/testing.md) guide to build the project and
  run its test suites.
- See [Writing Documentation](developer-guide/documentation.md) to learn how to document your code
  and regenerate this site.
- See [Linting & Formatting](developer-guide/linting.md) for the `format`/`docs-format`/`doc-lint`
  CMake targets that keep code and docs consistent.
- Browse the [API Reference](gecko/annotated.md) for the full class and namespace reference,
  generated directly from the Doxygen comments in the source code.
- Using Gecko from Python? See the [Python Bindings](user-guide/python.md) user guide.

## Project layout

Gecko is split into small, single-purpose CMake modules, each with its own `inc/` (public headers)
and `tst/` (Catch2 unit tests) directories:

| Module     | Depends on                                | Content                                                               |
| ---------- | ----------------------------------------- | --------------------------------------------------------------------- |
| `utils`    | —                                         | Strong ids, generic utilities (`gecko::StrongId`, ...)                |
| `math`     | `utils`                                   | `Point3d`, `Vector3d`, `BezierCurve`, `BezierSurface`, `BezierVolume` |
| `geom_itf` | `utils`, `math`                           | Concepts describing geometric entities                                |
| `geom`     | `utils`, `math`, `geom_itf`               | CAD model entities (`GeomVertex`, `GeomCurve`, ...)                   |
| `mesh`     | `utils`, `math`                           | `UnstructuredMesh`, `VariableRegistry`                                |
| `io`       | `utils`, `mesh`                           | Mesh structure file I/O (`GmshMeshReader`, `GmshMeshWriter`)          |
| `block`    | `utils`, `math`, `geom_itf`, `mesh`, `io` | CGAL-combinatorial-map-based quad/hex blocking (`Blocking`)           |

`gmds_core` also exists in the repository but is currently disabled in the root `CMakeLists.txt`
(its `add_subdirectory` call is commented out) and is not part of the build.

## Sandbox

`sandbox/` is a scratch executable (`gecko_sandbox`) for local development — a small, throwaway
`main.cpp` you're meant to edit freely while trying things out, built by default
(`-DGECKO_BUILD_SANDBOX=OFF` to skip it). It's not installed, not linked against by any other
module, and — unlike the modules above — not covered by the `format`/`doc-lint`/API reference
tooling: nothing in it is public API.

## Python bindings

`python/` builds a `gecko` Python extension module (pybind11), opt-in via `-DGECKO_BUILD_PYTHON=ON`
— currently an early scaffold, growing over time. See the [Python Bindings](user-guide/python.md)
guide for how to build it, use it, and run its test suite.
