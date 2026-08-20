# Linting & Formatting

Three independent, opt-in CMake targets keep the C++ code, the Markdown docs, and the Doxygen
comments themselves consistent. Each has a `-check`-style variant that only reports problems
(exit code, no changes — suitable for CI) and, where it makes sense, an apply variant that
rewrites files in place.

| Concern                               | Tool                                                                             | Apply           | Check only          |
| ------------------------------------- | -------------------------------------------------------------------------------- | --------------- | ------------------- |
| C++ code style                        | [clang-format](https://clang.llvm.org/docs/ClangFormat.html)                     | `format`        | `format-check`      |
| Markdown style (`docs/`, `README.md`) | [mdformat](https://mdformat.readthedocs.io/) + `mdformat-gfm` + `mdformat-admon` | `docs-format`   | `docs-format-check` |
| Doxygen comment completeness          | [Doxygen](https://www.doxygen.nl/) warnings (`Doxyfile.lint`)                    | — (no auto-fix) | `doc-lint`          |

Run everything check-only (no file changes) with the umbrella target:

```bash
cmake --build build --target lint
```

## Install the tools

```bash
brew install clang-format doxygen   # or: apt install clang-format doxygen
pip install -r requirements-docs.txt
```

Each target (controlled by `GECKO_BUILD_LINT`, `ON` by default) is only created when its
underlying tool is found on `PATH` at configure time; missing tools produce a `cmake` status
message rather than a hard failure, and `lint` itself only exists once all three are available.

## C++ formatting

```bash
cmake --build build --target format-check   # CI: fails if anything is unformatted
cmake --build build --target format         # rewrites files in place
```

The style is defined in `.clang-format` at the repository root (4-space indent, right-aligned
`&`/`*`, one-statement-per-line, ...) — run `clang-format -i --style=file <file>` on a single file
if you don't want to reformat everything. It covers the headers and tests of the library modules
(`utils`, `math`, `geom_itf`, `geom`, `mesh`, `io`, `block`), plus every source file under
`python/src/` and `biy/src/` — the whole first-party C++ tree.

## Markdown formatting

```bash
cmake --build build --target docs-format-check
cmake --build build --target docs-format
```

Configuration lives in `.mdformat.toml` at the repository root (keeps list numbers as typed,
preserves existing line wrapping). It reformats `docs/*.md` and `README.md`; MkDoxy's
auto-generated API Reference pages are not part of the repository and are never touched.

## Doxygen comment completeness

```bash
cmake --build build --target doc-lint
```

This runs Doxygen against `Doxyfile.lint` (at the repository root, `INPUT` restricted to the
library modules' `inc/` directories — `python/src/` and `biy/src/` are internal glue, not public
API, and are excluded here even though `format`/`format-check` above do cover them) with
`WARN_IF_UNDOCUMENTED`, `WARN_IF_DOC_ERROR` and `WARN_NO_PARAMDOC` enabled and
`WARN_AS_ERROR = FAIL_ON_WARNINGS`, so it fails the build on:

- undocumented public classes/structs/members/functions,
- `@param`/`@tparam` names that don't match the actual signature,
- other malformed `@`-commands.

It generates no HTML/site output (that's what `docs`/`docs-serve` are for, see
[Writing Documentation](documentation.md)) — it exists purely to catch the kind of comment/code
drift that's easy to introduce by hand, such as writing `@param AName` inside a `@brief` sentence
instead of using `@p AName` for an inline reference.
