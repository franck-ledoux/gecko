# -----------------------------------------------------------------------------
# Header-only library
# -----------------------------------------------------------------------------
add_library(gecko_block_lib INTERFACE)
# Alias that allows to write "target_link_libraries(my_exe PRIVATE Gecko::Block)"
add_library(Gecko::Block ALIAS gecko_block_lib)
# Directory where to find header files
target_include_directories(gecko_block_lib INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:include>
)

# List of all the header files to consider
target_sources(gecko_block_lib INTERFACE
        FILE_SET HEADERS
        BASE_DIRS inc
        FILES
        inc/gecko/block/Blocking.h
        inc/gecko/block/CellData.h
)

# Dependency to other Gecko modules
target_link_libraries(gecko_block_lib INTERFACE Gecko::Utils Gecko::Math Gecko::GeomItf Gecko::Mesh Gecko::Io)
# CGAL::Combinatorial_map (GPL-3.0-or-later, combinable with Gecko's AGPL-3.0 license under
# GPLv3/AGPLv3 section 13) is the backing data structure.
# BUILD_INTERFACE-only: cgal_dep/boost_dep are vendored, non-installed INTERFACE libraries, so they
# can't be part of gecko_block_lib's install(EXPORT ...) interface (see install() below).
target_link_libraries(gecko_block_lib INTERFACE $<BUILD_INTERFACE:cgal_dep> $<BUILD_INTERFACE:boost_dep>)
# Require C++23
target_compile_features(gecko_block_lib INTERFACE cxx_std_23)

# -----------------------------------------------------------------------------
# Test directory
# -----------------------------------------------------------------------------
if (GECKO_BUILD_TESTS)
    add_subdirectory(tst)
endif ()
# -----------------------------------------------------------------------------
# Exportation for external projects
# -----------------------------------------------------------------------------
# Allows another project to write: find_package(Gecko::Block REQUIRED)
install(TARGETS gecko_block_lib
        EXPORT GeckoBlockTargets
        FILE_SET HEADERS DESTINATION include
)

install(EXPORT GeckoBlockTargets
        FILE GeckoBlockTargets.cmake
        NAMESPACE Gecko::
        DESTINATION lib/cmake/GeckoBlock
)
