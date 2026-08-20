cmake_minimum_required(VERSION 3.20)
project(gecko VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#==============================================================================
# Operating System specificities (must be put at first)
#==============================================================================
if(APPLE)
    set(CMAKE_OSX_ARCHITECTURES arm64)
endif()
# -----------------------------------------------------------------------------
# Options du projet
# -----------------------------------------------------------------------------
# Default: ON, can be desactivated with -DGECKO_BUILD_TESTS=OFF
option(GECKO_BUILD_TESTS "Build the unit test of Gecko components" ON)
# Builds the sandbox/ scratch executable — a throwaway playground for local development, not
# installed and not part of the public API/docs. Default: ON, can be desactivated with
# -DGECKO_BUILD_SANDBOX=OFF
option(GECKO_BUILD_SANDBOX "Build the gecko_sandbox scratch executable" ON)
# Builds the python/ pybind11 module and its pytest suite. Default: OFF, since it needs a Python
# development environment (Python.h) that isn't guaranteed present; CI turns it ON explicitly.
option(GECKO_BUILD_PYTHON "Build the gecko_py Python bindings and their tests" OFF)
# Builds the biy/ interactive 3D blocking viewer. Default: OFF, since on top of a Python
# development environment it needs a GL/windowing stack (Polyscope's own dependency) that CI
# runners don't currently provide.
option(GECKO_BUILD_BIY "Build the biy interactive 3D blocking viewer" OFF)
# Option to activate/desactivate code coverage during unit testing
option(GECKO_CODE_COVERAGE "Perform code coverage with avail. test suites"               OFF)
option(GECKO_CODE_COVERAGE_HTML_REPORT "Activate html reporting during code covering"    OFF)

#==============================================================================
if (GECKO_CODE_COVERAGE_HTML_REPORT)
    set(GECKO_CODE_COVERAGE ON)
endif()
if (GECKO_CODE_COVERAGE)
    set(GECKO_TEST ON)
endif()

# =====================================================
# Embedded dependencies (Header-Only)
# =====================================================
# Boost
add_library(boost_dep INTERFACE)
target_include_directories(boost_dep INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/boost-1.87.0)

# CGAL
add_library(cgal_dep INTERFACE)
target_include_directories(cgal_dep INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/CGAL-5.6.2/include)

# Eigen3
add_library(eigen3_dep INTERFACE)
target_include_directories(eigen3_dep INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen-5.0.0)
add_library(Eigen3::Eigen ALIAS eigen3_dep)
# nlohmann_json
add_library(json_dep INTERFACE)
target_include_directories(json_dep INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/nlohmann_json/single_include)
# Note : allow to write #include <nlohmann/json.hpp> in the code


# =====================================================
# Embedded dependencies (Compiled)
# =====================================================

# Catch2 v3
set(CATCH_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(CATCH_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(external/Catch2)

# pybind11 (only pulled in — and only requires a Python dev environment — when actually needed).
# biy needs it too, for its embedded interpreter.
if (GECKO_BUILD_PYTHON OR GECKO_BUILD_BIY)
    add_subdirectory(external/pybind11)
endif ()

# Polyscope (biy's 3D viewer). Vendored with its own dependencies (GLFW, ImGui, glm, ...), but
# still needs the platform's GL/windowing headers, hence its own opt-in option.
if (GECKO_BUILD_BIY)
    add_subdirectory(external/polyscope)
endif ()


# =====================================================
# Test configuration
# =====================================================
if (GECKO_BUILD_TESTS)
    enable_testing() # Active CTest

    set(TEST_SAMPLES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/test_data")
    configure_file(
            unit_test_config.h.in
            ${CMAKE_CURRENT_BINARY_DIR}/unit_test_config.h
            @ONLY
    )
endif ()


#==============================================================================
# TESTING AND CODE COVERAGE
#
# Up to now, some components are not reliable enough to be added in the full
# code coverage process. It must be improved in the next few months
#==============================================================================
if (GECKO_CODE_COVERAGE)
    message(STATUS "============================================")
    message(STATUS "============= COVERAGE MODE  ===============")
    message(STATUS "============================================")
    message(STATUS "CXX Compiler          = ${CMAKE_CXX_COMPILER}")
    message(STATUS "CXX Compiler Version  = ${CMAKE_CXX_COMPILER_VERSION}")
    message(STATUS "CMAKE_CXX_FLAGS       = ${CMAKE_CXX_FLAGS}")
    message(STATUS "CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")
    message(STATUS "============================================")
    message(STATUS "============ COVERED COMPONENTS ============")
    message(STATUS "============================================")
    foreach(c IN LISTS GMDS_COVER_COMPONENTS)
        set(ENABLE_${c} ON)
        message(STATUS "COMPONENT  : ${c}")
    endforeach()
    message(STATUS "============================================")
endif(GECKO_CODE_COVERAGE)
#==============================================================================
#voir a ajouter une target coverage avec ces options liées à la target
# va m'éviter les variables locales
# target vide si pas besoin ou pas
# forcer mode non optimise plutôt
# Attention ne pas livrer avec code coverage!!
# Faire un build type Coverage, comme Debug ou Release
if (GECKO_CODE_COVERAGE)
    set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)
    include(CodeCoverage)
    set(CMAKE_CXX_FLAGS "-g -O0 --coverage -fprofile-arcs -ftest-coverage")
    set(CMAKE_C_FLAGS "-g -O0 --coverage -fprofile-arcs -ftest-coverage")
    set(COVERAGE_LCOV_EXCLUDES
            '/usr/*'
            '${CMAKE_SOURCE_DIR}/docs/*'
            '${CMAKE_SOURCE_DIR}/test_data/*'
            '${CMAKE_SOURCE_DIR}/external/*'
            '${CMAKE_SOURCE_DIR}/*/tst/*'
    )

    #==============================================================================
    # Data directory if required for some tests
    set(DATA_DIR ${CMAKE_BINARY_DIR}/test_data)
    # Directory that contains test suites

    if(GECKO_CODE_COVERAGE_HTML_REPORT)
        SETUP_TARGET_FOR_COVERAGE_LCOV(
                NAME code_cover_gecko
                EXECUTABLE ctest
                ${DATA_DIR} .
        )
    else(GECKO_CODE_COVERAGE_HTML_REPORT)
        SETUP_TARGET_FOR_COVERAGE_LCOV_SHORT(
                NAME code_cover_gecko
                EXECUTABLE ctest
                ${DATA_DIR}
        )
    endif(GECKO_CODE_COVERAGE_HTML_REPORT)
    #==============================================================================
    #voir configure_file qui garde le lien de mise à jour, dépendances conservées
    file(COPY
            ${CMAKE_SOURCE_DIR}/test_data/
            DESTINATION ${CMAKE_BINARY_DIR}/test_data/)

endif(GECKO_CODE_COVERAGE)


#==============================================================================
# DEVELOPER DOCUMENTATION (Doxygen comments -> MkDocs + Material + MkDoxy)
#==============================================================================
option(GECKO_BUILD_DOCS "Enable the 'docs'/'docs-serve' targets to build the MkDocs developer documentation" ON)
if (GECKO_BUILD_DOCS)
    find_program(DOXYGEN_EXECUTABLE NAMES doxygen)
    find_program(MKDOCS_EXECUTABLE NAMES mkdocs)

    if (DOXYGEN_EXECUTABLE AND MKDOCS_EXECUTABLE)
        # Builds the static site into ${CMAKE_SOURCE_DIR}/site
        add_custom_target(docs
                COMMAND ${MKDOCS_EXECUTABLE} build --clean
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                COMMENT "Generating developer documentation with MkDocs (see site/index.html)"
                VERBATIM
        )
        # Serves the site locally with live-reload on http://127.0.0.1:8000
        add_custom_target(docs-serve
                COMMAND ${MKDOCS_EXECUTABLE} serve
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                COMMENT "Serving developer documentation on http://127.0.0.1:8000"
                USES_TERMINAL
                VERBATIM
        )
    else ()
        message(STATUS "GECKO_BUILD_DOCS: doxygen and/or mkdocs not found on PATH, 'docs'/'docs-serve' targets disabled")
        message(STATUS "  -> install them, e.g.: brew install doxygen && pip install -r requirements-docs.txt")
    endif ()
endif (GECKO_BUILD_DOCS)


#==============================================================================
# LINTING (clang-format for C++, mdformat for Markdown, Doxygen for doc completeness)
#==============================================================================
option(GECKO_BUILD_LINT "Enable the 'format'/'format-check'/'docs-format'/'docs-format-check'/'doc-lint'/'lint' targets" ON)
if (GECKO_BUILD_LINT)
    set(GECKO_LINT_MODULES utils math geom_itf geom mesh io block python biy)

    file(GLOB_RECURSE GECKO_CXX_SOURCES CONFIGURE_DEPENDS
            ${CMAKE_SOURCE_DIR}/utils/inc/*.h ${CMAKE_SOURCE_DIR}/utils/tst/*.cpp
            ${CMAKE_SOURCE_DIR}/math/inc/*.h ${CMAKE_SOURCE_DIR}/math/tst/*.cpp
            ${CMAKE_SOURCE_DIR}/geom_itf/inc/*.h ${CMAKE_SOURCE_DIR}/geom_itf/tst/*.cpp
            ${CMAKE_SOURCE_DIR}/geom/inc/*.h ${CMAKE_SOURCE_DIR}/geom/tst/*.cpp
            ${CMAKE_SOURCE_DIR}/mesh/inc/*.h ${CMAKE_SOURCE_DIR}/mesh/tst/*.cpp
            ${CMAKE_SOURCE_DIR}/io/inc/*.h ${CMAKE_SOURCE_DIR}/io/tst/*.cpp
            ${CMAKE_SOURCE_DIR}/block/inc/*.h ${CMAKE_SOURCE_DIR}/block/tst/*.cpp
            ${CMAKE_SOURCE_DIR}/python/src/*.h ${CMAKE_SOURCE_DIR}/python/src/*.cpp
            ${CMAKE_SOURCE_DIR}/biy/src/*.h ${CMAKE_SOURCE_DIR}/biy/src/*.cpp
    )
    file(GLOB_RECURSE GECKO_MARKDOWN_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/docs/*.md)
    list(APPEND GECKO_MARKDOWN_SOURCES ${CMAKE_SOURCE_DIR}/README.md)

    # --- C++ formatting (clang-format), scoped to the modules above ---
    find_program(CLANG_FORMAT_EXECUTABLE NAMES clang-format)
    if (CLANG_FORMAT_EXECUTABLE)
        add_custom_target(format
                COMMAND ${CLANG_FORMAT_EXECUTABLE} -i --style=file ${GECKO_CXX_SOURCES}
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                COMMENT "Formatting C++ sources with clang-format (in place)"
                VERBATIM
        )
        add_custom_target(format-check
                COMMAND ${CLANG_FORMAT_EXECUTABLE} --dry-run --Werror --style=file ${GECKO_CXX_SOURCES}
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                COMMENT "Checking C++ formatting with clang-format (no changes made)"
                VERBATIM
        )
    else ()
        message(STATUS "GECKO_BUILD_LINT: clang-format not found on PATH, 'format'/'format-check' targets disabled")
        message(STATUS "  -> install it, e.g.: brew install clang-format")
    endif ()

    # --- Markdown formatting (mdformat), for docs/ and README.md ---
    find_program(MDFORMAT_EXECUTABLE NAMES mdformat)
    if (MDFORMAT_EXECUTABLE)
        add_custom_target(docs-format
                COMMAND ${MDFORMAT_EXECUTABLE} ${GECKO_MARKDOWN_SOURCES}
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                COMMENT "Formatting Markdown docs with mdformat (in place)"
                VERBATIM
        )
        add_custom_target(docs-format-check
                COMMAND ${MDFORMAT_EXECUTABLE} --check ${GECKO_MARKDOWN_SOURCES}
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                COMMENT "Checking Markdown formatting with mdformat (no changes made)"
                VERBATIM
        )
    else ()
        message(STATUS "GECKO_BUILD_LINT: mdformat not found on PATH, 'docs-format'/'docs-format-check' targets disabled")
        message(STATUS "  -> install it, e.g.: pip install -r requirements-docs.txt")
    endif ()

    # --- Doxygen comment completeness/consistency check (undocumented members, @param mismatches, ...) ---
    # (independent of GECKO_BUILD_DOCS, hence its own find_program call)
    find_program(DOXYGEN_EXECUTABLE NAMES doxygen)
    if (DOXYGEN_EXECUTABLE)
        add_custom_target(doc-lint
                COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile.lint
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                COMMENT "Checking Doxygen comment completeness (Doxyfile.lint)"
                VERBATIM
        )
    else ()
        message(STATUS "GECKO_BUILD_LINT: doxygen not found on PATH, 'doc-lint' target disabled")
    endif ()

    # --- Umbrella check-only target, suitable for CI ---
    if (CLANG_FORMAT_EXECUTABLE AND MDFORMAT_EXECUTABLE AND DOXYGEN_EXECUTABLE)
        add_custom_target(lint DEPENDS format-check docs-format-check doc-lint)
    endif ()
endif (GECKO_BUILD_LINT)


# =====================================================
# Modules du projet
# =====================================================
add_subdirectory(utils)
add_subdirectory(math)
add_subdirectory(geom_itf)
add_subdirectory(mesh)
add_subdirectory(io)
add_subdirectory(geom)
add_subdirectory(block)

if (GECKO_BUILD_SANDBOX)
    add_subdirectory(sandbox)
endif ()

if (GECKO_BUILD_PYTHON)
    add_subdirectory(python)
endif ()

if (GECKO_BUILD_BIY)
    add_subdirectory(biy)
endif ()