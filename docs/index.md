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
- Want to see (and build) a blocking interactively? See the [biy viewer](user-guide/biy.md) guide.

## Project layout

Gecko's library code is split into small, single-purpose CMake modules, each with its own `inc/`
(public headers) and `tst/` (Catch2 unit tests) directories, built in this dependency order:

| Module     | Depends on                                | Content                                                               |
| ---------- | ----------------------------------------- | --------------------------------------------------------------------- |
| `utils`    | —                                         | Strong ids, generic utilities (`gecko::StrongId`, ...)                |
| `math`     | `utils`                                   | `Point3d`, `Vector3d`, `BezierCurve`, `BezierSurface`, `BezierVolume` |
| `geom_itf` | `utils`, `math`                           | Concepts describing geometric entities                                |
| `mesh`     | `utils`, `math`                           | `UnstructuredMesh`, `VariableRegistry`                                |
| `io`       | `utils`, `mesh`                           | Mesh structure file I/O (`GmshMeshReader`, `GmshMeshWriter`)          |
| `geom`     | `utils`, `math`, `geom_itf`, `mesh`, `io` | CAD model entities (`FacetedGeometry`, `FacetedVertex`, ...)          |
| `block`    | `utils`, `math`, `geom_itf`, `mesh`, `io` | CGAL-combinatorial-map-based quad/hex blocking (`Blocking`)           |

On top of those, two more directories build against the modules above but aren't themselves
libraries other code links against:

| Directory | Depends on                                                               | Content                                                                             |
| --------- | ------------------------------------------------------------------------ | ----------------------------------------------------------------------------------- |
| `python`  | `Gecko::Block`, `Gecko::Geom`                                            | The `gecko` pybind11 extension module — see [Python Bindings](user-guide/python.md) |
| `biy`     | `Gecko::Block`, `Gecko::Geom`, `Gecko::Io`, Polyscope, `pybind11::embed` | The interactive 3D viewer executable — see [the biy guide](user-guide/biy.md)       |

## Sandbox

`sandbox/` is a scratch executable (`gecko_sandbox`) for local development — a small, throwaway
`main.cpp` you're meant to edit freely while trying things out, built by default
(`-DGECKO_BUILD_SANDBOX=OFF` to skip it). It's not installed, not linked against by any other
module, and — unlike the modules above — not covered by the `format`/`doc-lint`/API reference
tooling: nothing in it is public API.

## Python bindings

`python/` builds a `gecko` Python extension module (pybind11), opt-in via `-DGECKO_BUILD_PYTHON=ON`
— façade managers (`GeomModel`, `Blocking`) over the C++ API, growing over time. See the
[Python Bindings](user-guide/python.md) guide for how to build it, use it, and run its test suite.

## biy — interactive 3D viewer

`biy/` builds `biy` ("Block It Yourself"), an interactive [Polyscope](https://polyscope.run)-based
3D viewer for building a `Blocking` against a `GeomModel`, opt-in via `-DGECKO_BUILD_BIY=ON`. It
reuses the same `GeomModel`/`Blocking` façade as the Python bindings above — both an on-screen panel
and an embedded Python console drive the exact same live objects. See
[the biy guide](user-guide/biy.md) for how to build and use it.

!!! warning
    Unlike every other target here, `biy` is not built or tested in CI: the hosted runners provide
    no GL/windowing stack for Polyscope. See [the biy guide](user-guide/biy.md#building) for details.
