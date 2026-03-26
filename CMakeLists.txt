cmake_minimum_required(VERSION 3.20)

project(gecko VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

include(FetchContent)

# =====================================================
# Bundled dependencies (physically present in subprojects/)
# =====================================================

# Boost 1.87.0 (bundled, header-only)
add_library(boost_dep INTERFACE)
target_include_directories(boost_dep INTERFACE
    ${CMAKE_SOURCE_DIR}/subprojects/boost-1.87.0
)

# CGAL 5.6.2 (bundled, header-only)
add_library(cgal_dep INTERFACE)
target_include_directories(cgal_dep INTERFACE
    ${CMAKE_SOURCE_DIR}/subprojects/CGAL-5.6.2/include
)

# predicates (bundled, compiled)
add_subdirectory(subprojects/predicates)

# =====================================================
# Downloaded dependencies (Meson wrap equivalents)
# =====================================================

# Eigen3 3.4.0
FetchContent_Declare(
    eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2
    URL_HASH SHA256=b4c198460eba6f28d34894e3a5710998818515104d6e74e5cc331ce31e46e626
)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(eigen3)

# nlohmann_json 3.11.3
FetchContent_Declare(
    nlohmann_json
    URL https://github.com/nlohmann/json/archive/refs/tags/v3.11.3.tar.gz
    URL_HASH SHA256=0d8ef5af7f9794e3263480193c491549b2ba6cc74bb018906202ada498a79406
)
set(JSON_BuildTests OFF CACHE BOOL "" FORCE)
set(JSON_Install OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(nlohmann_json)

# Catch2 3.8.0
FetchContent_Declare(
    Catch2
    URL https://github.com/catchorg/Catch2/archive/v3.8.0.tar.gz
    URL_HASH SHA256=1ab2de20460d4641553addfdfe6acd4109d871d5531f8f519a52ea4926303087
)
FetchContent_MakeAvailable(Catch2)

# =====================================================
# Unit test config header (equivalent to Meson configure_file)
# =====================================================
set(TEST_SAMPLES_DIR "${CMAKE_SOURCE_DIR}/test_data")
configure_file(
    ${CMAKE_SOURCE_DIR}/unit_test_config.h.in
    ${CMAKE_BINARY_DIR}/unit_test_config.h
    @ONLY
)

# =====================================================
# Project modules
# =====================================================
add_subdirectory(cblock)

# =====================================================
# Tests
# =====================================================
enable_testing()
