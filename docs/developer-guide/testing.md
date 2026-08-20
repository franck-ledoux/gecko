# Building & Testing

## Prerequisites

- CMake >= 3.20
- A C++20/23 compiler (GCC, Clang, MSVC/AppleClang)
- [Catch2](https://github.com/catchorg/Catch2) v3 — bundled under `external/Catch2`, built automatically
- For code coverage reports: [`lcov`](https://github.com/linux-test-project/lcov) and `genhtml`
  (installed together, e.g. `brew install lcov` on macOS, `apt install lcov` on Debian/Ubuntu)

## Building the project

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
```

`GECKO_BUILD_TESTS` is `ON` by default and registers every module's test executable with CTest.
Pass `-DGECKO_BUILD_TESTS=OFF` to skip building the test suites.

## Running the tests

Once built, run every test suite with CTest from the build directory:

```bash
cd build
ctest --output-on-failure
```

Run a single module's suite directly, or filter CTest by name:

```bash
./utils/tst/test_utils
ctest -R test_math -V
```

Each module registers exactly one CTest entry per test executable (see `<module>/tst/CMakeLists.txt`),
so `ctest -R <name>` matches the executable name (`test_utils`, `test_math`, `test_geometry`, `test_mesh`, ...).

## Writing a new test

Tests use [Catch2](https://github.com/catchorg/Catch2) v3. To add a test to an existing module:

1. Add a `.cpp` file under `<module>/tst/`, `#include <catch2/catch_test_macros.hpp>` (or the
   relevant Catch2 header for the matchers/generators you need), and write `TEST_CASE`s.

2. Register the new file in `<module>/tst/CMakeLists.txt`, e.g.:

   ```cmake
   add_executable(test_math
           bezier_curve_tests.cpp
           PointVector_tests.cpp
           my_new_tests.cpp)          # <- add your file here
   ```

   If the module doesn't have a test executable yet, follow the pattern used by the existing
   modules (`utils/tst/CMakeLists.txt` is the simplest example):

   ```cmake
   add_executable(test_utils utils_tests.cpp)
   target_link_libraries(test_utils PRIVATE
       Gecko::Utils
       Catch2::Catch2
       Catch2::Catch2WithMain
   )
   target_include_directories(test_utils PRIVATE ${CMAKE_BINARY_DIR})
   add_test(NAME test_utils COMMAND test_utils)
   ```

3. Re-run CMake configure (so the new executable/test is picked up) and build:

   ```bash
   cmake -S . -B build
   cmake --build build --target test_math -j
   ```

## Code coverage (HTML report)

Coverage is driven by two CMake options, both `OFF` by default:

- `GECKO_CODE_COVERAGE` — compiles with `--coverage` and enables the `code_cover_gecko` target
  (`lcov`'s compact text summary).
- `GECKO_CODE_COVERAGE_HTML_REPORT` — same as above, plus renders an HTML report with `genhtml`.

Use a **separate build directory** for coverage builds, since they change compiler/linker flags:

```bash
cmake -S . -B build-coverage -DCMAKE_BUILD_TYPE=Debug -DGECKO_CODE_COVERAGE_HTML_REPORT=ON
cmake --build build-coverage -j
cmake --build build-coverage --target code_cover_gecko
```

The `code_cover_gecko` target resets counters, runs `ctest`, then produces:

- `build-coverage/code_cover_gecko.info` — the raw lcov trace file.
- `build-coverage/code_cover_gecko/index.html` — the HTML report; open it in a browser.

!!! note
    Coverage runs with `-DGECKO_BUILD_PYTHON=ON` (see `.github/workflows/coverage.yml`), so the
    library modules (`utils`, `math`, `geom_itf`, `mesh`, `geom`, `io`, `block`) and the `python/`
    bindings are all covered — running the Python test suite exercises the same instrumented C++
    build as the C++ suites do (see [Python Bindings](../user-guide/python.md#code-coverage)).
    `biy/` is excluded, since it isn't built in CI at all (see
    [Project layout](../index.md#project-layout)).
