# Python Bindings

Gecko ships a `gecko` Python extension module built with [pybind11](https://pybind11.readthedocs.io/).
It's currently an early scaffold — enough to prove the build/test/coverage/documentation pipeline
works end-to-end — rather than a full binding of the C++ API; the surface area will grow from here.

!!! note
    This page, and the example on it, is exercised by [`python/tst/test_hello.py`](https://github.com/franck-ledoux/gecko/blob/main/python/tst/test_hello.py),
    which CI runs on every push. If this page and the actual behavior ever drift apart, that test
    (and therefore CI) fails — so what's documented here is always what the code actually does.

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
