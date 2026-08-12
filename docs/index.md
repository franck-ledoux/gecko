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

## Project layout

Gecko is split into small, single-purpose CMake modules, each with its own `inc/` (public headers)
and `tst/` (Catch2 unit tests) directories:

| Module     | Depends on                  | Content                                                |
| ---------- | --------------------------- | ------------------------------------------------------ |
| `utils`    | —                           | Strong ids, generic utilities (`gecko::StrongId`, ...) |
| `math`     | `utils`                     | `Point3d`, `Vector3d`, `BezierCurve`                   |
| `geom_itf` | `utils`, `math`             | Concepts describing geometric entities                 |
| `geom`     | `utils`, `math`, `geom_itf` | CAD model entities (`GeomVertex`, `GeomCurve`, ...)    |
| `mesh`     | `utils`, `math`             | `UnstructuredMesh`, `VariableRegistry`                 |

`block` and `gmds_core` also exist in the repository but are currently disabled in the root
`CMakeLists.txt` (`add_subdirectory` calls are commented out) and are not part of the build yet.
