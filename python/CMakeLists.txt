# -----------------------------------------------------------------------------
# gecko: pybind11 Python bindings over Gecko::App's façades, which are what actually wrap the
# kernel in primitive types; this module is only the pybind11 glue. See docs/user-guide/python.md.
# Behind GECKO_BUILD_PYTHON (root CMakeLists.txt).
# -----------------------------------------------------------------------------
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

pybind11_add_module(gecko
        src/bindings.cpp
        src/PyGeckoBindings.cpp
)
target_link_libraries(gecko PRIVATE Gecko::App)
# Matches Gecko::Block/Gecko::Geom's own required C++ standard.
target_compile_features(gecko PRIVATE cxx_std_23)

if (GECKO_BUILD_TESTS)
    add_test(
            NAME test_python
            COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/tst -v
    )
    # PYTHONPATH points at the just-built extension's own directory, so `import gecko` finds it
    # without needing a pip install step. GECKO_TEST_DATA_DIR lets Python tests locate the same
    # checked-in .msh fixtures the C++ suites use via unit_test_config.h's TEST_SAMPLES_DIR.
    set_tests_properties(test_python PROPERTIES
            ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:gecko>;GECKO_TEST_DATA_DIR=${CMAKE_SOURCE_DIR}/test_data"
    )
endif ()
