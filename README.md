# GECKO

GECKO stands for hiGh ordEr bloCK tOolkit. It provides data structures and algorithms to easily modify a block structure used for discretizing a geometrical domain with a hexahedral mesh.

## Prerequisites

- CMake >= 3.20
- A C++17 compiler (GCC, Clang, MSVC)
- An internet connection (CMake will automatically download Eigen3, nlohmann_json and Catch2)

The following dependencies are bundled in the repository under `subprojects/`:
- Boost 1.87.0 (header-only)
- CGAL 5.6.2 (header-only)
- predicates (compiled)

## Build

```sh
mkdir build && cd build
cmake ..
make
```

## Run unit tests

Tests are located in `cblock/tst/`. After building, run them with:

```sh
cd build
ctest -R test_cblock -V
```

To rebuild and run in one step:

```sh
cd build
make test_cblock && ctest -R test_cblock -V
```
