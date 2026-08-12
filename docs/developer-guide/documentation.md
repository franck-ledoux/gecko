# Writing Documentation

This site is generated from the Doxygen-style comments written directly in the public headers
(`<module>/inc/gecko/...`). There is nothing to write by hand outside of the source code and the
handful of Markdown pages under `docs/` (this page included) — the [API Reference](../gecko/annotated.md)
is produced automatically from your comments every time the site is rebuilt.

## Comment style

Gecko follows the classic Doxygen `@`-command style (Javadoc-like), illustrated in
`math/inc/gecko/math/BezierCurve.h` — use that file as the reference example when documenting new
code. In short:

- Every class/struct gets a `@class`/`@brief` block above its declaration, plus `@tparam` for each
  template parameter.
- Every public method gets a `/** ... */` block with:
  - `@brief` — a one-line summary of what it does.
  - `@tparam` — one per template parameter.
  - `@param` — one per function parameter.
  - `@return` — what the function returns (omit for `void`).
  - `@pre` — preconditions the caller must satisfy, if any.
  - `@throw` — exceptions that can be thrown, if any.
  - `@see` — cross-references to related functions, if useful.
- One-liners (simple accessors, thin operators) can use the compact `/** @brief ... */` form on a
  single line instead of a full multi-line block.

```cpp
/**
 * @brief Evaluates the point on the curve at parameter t using De Casteljau's algorithm.
 * @param AC curvilinear coordinate in the range [0.0, 1.0].
 * @return The evaluated point of type PointT at parameter t.
 * @pre Parameter AC must satisfy 0.0 <= t <= 1.0 (enforced via assert in Debug builds).
 */
TPointT value(double AC) const { /* ... */ }
```

Keep comments focused on the **why/contract** (preconditions, ownership, complexity, units) rather
than restating what a well-named function obviously does.

## Toolchain

The site is built with:

- [Doxygen](https://www.doxygen.nl/) — parses the C++ headers and the `@`-comments into an XML model.
- [MkDocs](https://www.mkdocs.org/) + [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/) — renders the site.
- [MkDoxy](https://mkdoxy.kubaandrysek.cz/) — the MkDocs plugin that turns the Doxygen XML into the
  Material-themed API Reference pages, and drives Doxygen itself (no separate `Doxyfile` to maintain
  — its configuration lives in `mkdocs.yml` at the repository root, under `plugins.mkdoxy`).

### Install the toolchain

```bash
brew install doxygen        # or: apt install doxygen / dnf install doxygen
pip install -r requirements-docs.txt
```

`requirements-docs.txt` pins `mkdocs`, `mkdocs-material` and `mkdoxy`.

### Generate the site

A CMake target wraps `mkdocs build` so documentation generation fits the same workflow as building
and testing:

```bash
cmake -S . -B build
cmake --build build --target docs
```

This writes a static site to `site/` at the repository root; open `site/index.html` in a browser.

### Live preview while writing

`mkdocs serve` rebuilds on every save (Markdown pages **and** header comments) and serves the site
locally:

```bash
cmake --build build --target docs-serve
# or directly: mkdocs serve
```

Then browse to <http://127.0.0.1:8000>.

!!! note
    The `docs`/`docs-serve` CMake targets are only created when both `doxygen` and `mkdocs` are
    found on the `PATH` at configure time (controlled by the `GECKO_BUILD_DOCS` option, `ON` by
    default). If they're missing, `cmake` prints a status message and skips the targets — install
    the toolchain above and re-run `cmake` to enable them.

## Adding a new documentation page

Hand-written pages (like this one) live under `docs/` as Markdown files, and must be listed in the
`nav:` section of `mkdocs.yml` to appear in the site's navigation. The API Reference pages
(`gecko/annotated.md`, `gecko/classes.md`, ...) are generated automatically by MkDoxy and should
not be created or edited by hand.
