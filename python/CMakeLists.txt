# -----------------------------------------------------------------------------
# gecko: pybind11 Python bindings — currently a minimal "hello world" scaffold proving the infra
# (build, pytest via CTest, code coverage of both the bindings and the C++ library, docs) works
# end-to-end; see docs/user-guide/python.md. Behind GECKO_BUILD_PYTHON (root CMakeLists.txt).
# -----------------------------------------------------------------------------
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

pybind11_add_module(gecko src/bindings.cpp)
target_compile_features(gecko PRIVATE cxx_std_20)

if (GECKO_BUILD_TESTS)
    add_test(
            NAME test_python
            COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/tst -v
    )
    # PYTHONPATH points at the just-built extension's own directory, so `import gecko` finds it
    # without needing a pip install step.
    set_tests_properties(test_python PROPERTIES
            ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:gecko>"
    )
endif ()
