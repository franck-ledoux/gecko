# -----------------------------------------------------------------------------
# Header-only library
# -----------------------------------------------------------------------------
add_library(gecko_geom_lib INTERFACE)
# Alias that allows to write "target_link_libraries(my_exe PRIVATE Gecko::Geom)"
add_library(Gecko::Geom ALIAS gecko_geom_lib)
# Directory where to find header files
target_include_directories(gecko_geom_lib INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:include>
)

# List of all the header files to consider
target_sources(gecko_geom_lib INTERFACE
        FILE_SET HEADERS
        BASE_DIRS inc
        FILES
        inc/gecko/geom/GeomEntityRefVariant.h
        inc/gecko/geom/GeomEntities.h
        inc/gecko/geom/Geometry.h
)

# Dependency to other Gecko modules
target_link_libraries(gecko_geom_lib INTERFACE Gecko::Utils Gecko::Math Gecko::GeomItf)
# Require C++23
target_compile_features(gecko_geom_lib INTERFACE cxx_std_23)
# -----------------------------------------------------------------------------
# Test directory
# -----------------------------------------------------------------------------
if (GECKO_BUILD_TESTS)
    add_subdirectory(tst)
endif ()
# -----------------------------------------------------------------------------
# Exportation for external projects
# -----------------------------------------------------------------------------
# Allows another project to write: find_package(Gecko::Geom REQUIRED)
install(TARGETS gecko_geom_lib
        EXPORT GeckoGeomTargets
        FILE_SET HEADERS DESTINATION include
)

install(EXPORT GeckoGeomTargets
        FILE GeckoGeomTargets.cmake
        NAMESPACE Gecko::
        DESTINATION lib/cmake/GeckoGeom
)
