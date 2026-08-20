# GECKO

[![Ubuntu](https://img.shields.io/github/actions/workflow/status/franck-ledoux/gecko/ci.yml?branch=main&event=push&job=build-test%20%28ubuntu-latest%29&label=Ubuntu)](https://github.com/franck-ledoux/gecko/actions/workflows/ci.yml)
[![macOS](https://img.shields.io/github/actions/workflow/status/franck-ledoux/gecko/ci.yml?branch=main&event=push&job=build-test%20%28macos-latest%29&label=macOS)](https://github.com/franck-ledoux/gecko/actions/workflows/ci.yml)
[![Windows](https://img.shields.io/github/actions/workflow/status/franck-ledoux/gecko/ci.yml?branch=main&event=push&job=build-test%20%28windows-latest%29&label=Windows)](https://github.com/franck-ledoux/gecko/actions/workflows/ci.yml)
[![Coverage](https://github.com/franck-ledoux/gecko/actions/workflows/coverage.yml/badge.svg?branch=main)](https://github.com/franck-ledoux/gecko/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/franck-ledoux/gecko/graph/badge.svg?token=MA8EPY7I94)](https://codecov.io/gh/franck-ledoux/gecko)

GECKO stands for hi**G**h ord**E**r blo**C**K t**O**olkit. It provides data structures and
algorithms to easily modify a block structure used for discretizing a geometrical domain with a
hexahedral mesh.

## Quick start

Requires CMake >= 3.20 and a C++20/23 compiler (GCC, Clang, MSVC/AppleClang).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
cd build && ctest --output-on-failure
```

## biy — Block It Yourself

`biy` is Gecko's interactive 3D viewer for building a block structure against a geometric model.
It loads a Gmsh `.msh` file, shows it, and lets you create a bounding box, drag its corners onto the
geometry, mesh it and classify it — from the on-screen panel, with the mouse, or by typing Python
at a console that drives the very same live objects.

![The biy viewer: a hex block fitted inside a translucent cylinder, its corners colored by what
they are classified on](docs/assets/biy.png)

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DGECKO_BUILD_BIY=ON
cmake --build build --target biy
./build/biy/biy test_data/cylinder.msh
```

It is opt-in (`GECKO_BUILD_BIY=OFF` by default) because, on top of a Python development
environment, it needs a GL/windowing stack. **CI does not build or test `biy`**, for the same
reason — the hosted runners provide no such stack. Its supporting library code (`Blocking`, the
Python facades) *is* covered by the usual suites; the viewer itself is verified locally.

See [the biy guide](docs/user-guide/biy.md) for the mouse modes, the corner color code and the
`biy_config.json` settings.

## Documentation

The full developer documentation — project layout, running/writing tests, code coverage,
writing/generating this documentation itself, and linting — lives under [`docs/`](docs/index.md)
and is kept up to date there rather than duplicated here:

- [Overview & project layout](docs/index.md)
- [biy — the interactive 3D viewer](docs/user-guide/biy.md)
- [Python bindings](docs/user-guide/python.md)
- [Building & Testing](docs/developer-guide/testing.md)
- [Writing Documentation](docs/developer-guide/documentation.md)
- [Linting & Formatting](docs/developer-guide/linting.md)

For the full experience (searchable, with the API reference generated from the source code
comments), browse the published site — **<https://franck-ledoux.github.io/gecko/>**, rebuilt from
`main` on every push — instead of reading the raw Markdown. To build and browse it locally instead:

```bash
cmake --build build --target docs-serve   # http://127.0.0.1:8000, live-reloading
# or: cmake --build build --target docs   # static site in site/index.html
```
