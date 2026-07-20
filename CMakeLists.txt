cmake_minimum_required(VERSION 3.20)
project(gecko VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

# predicates
add_subdirectory(external/predicates)

# Catch2 v3 (Remplace FetchContent)
set(CATCH_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(CATCH_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(external/catch2)


# =====================================================
# Test configuration
# =====================================================
enable_testing()

set(TEST_SAMPLES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/test_data")
configure_file(
        unit_test_config.h.in
        ${CMAKE_CURRENT_BINARY_DIR}/unit_test_config.h
        @ONLY
)


# =====================================================
# Modules du projet
# =====================================================
add_subdirectory(gmds_core)
add_subdirectory(block)
