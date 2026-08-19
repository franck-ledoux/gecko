# GECKO

**Build & Test:** Ubuntu
[![Ubuntu](https://github.com/franck-ledoux/gecko/actions/workflows/ci.yml/badge.svg?branch=main&job=build-test%20%28ubuntu-latest%29)](https://github.com/franck-ledoux/gecko/actions/workflows/ci.yml)
macOS
[![macOS](https://github.com/franck-ledoux/gecko/actions/workflows/ci.yml/badge.svg?branch=main&job=build-test%20%28macos-latest%29)](https://github.com/franck-ledoux/gecko/actions/workflows/ci.yml)
Windows
[![Windows](https://github.com/franck-ledoux/gecko/actions/workflows/ci.yml/badge.svg?branch=main&job=build-test%20%28windows-latest%29)](https://github.com/franck-ledoux/gecko/actions/workflows/ci.yml)
<br>
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

## Documentation

The full developer documentation — project layout, running/writing tests, code coverage,
writing/generating this documentation itself, and linting — lives under [`docs/`](docs/index.md)
and is kept up to date there rather than duplicated here:

- [Overview & project layout](docs/index.md)
- [Building & Testing](docs/developer-guide/testing.md)
- [Writing Documentation](docs/developer-guide/documentation.md)
- [Linting & Formatting](docs/developer-guide/linting.md)

For the full experience (searchable, with the API reference generated from the source code
comments), build and browse the site instead of reading the raw Markdown:

```bash
cmake --build build --target docs-serve   # http://127.0.0.1:8000, live-reloading
# or: cmake --build build --target docs   # static site in site/index.html
```
