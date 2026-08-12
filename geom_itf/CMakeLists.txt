# -----------------------------------------------------------------------------
# Header-only library
# -----------------------------------------------------------------------------
add_library(gecko_geom_itf_lib INTERFACE)
# Alias that allows to write "target_link_libraries(my_exe PRIVATE Gecko::GeomItf)"
add_library(Gecko::GeomItf ALIAS gecko_geom_itf_lib)
# Directory where to find header files
target_include_directories(gecko_geom_itf_lib INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:include>
)

# List of all the header files to consider
target_sources(gecko_geom_itf_lib INTERFACE
        FILE_SET HEADERS
        BASE_DIRS inc
        FILES
        inc/gecko/geom_itf/Concepts.h
)

# Dependency to other Gecko modules
target_link_libraries(gecko_geom_itf_lib INTERFACE Gecko::Utils Gecko::Math)
# Require C++23
target_compile_features(gecko_geom_itf_lib INTERFACE cxx_std_23)

# -----------------------------------------------------------------------------
# Exportation for external projects
# -----------------------------------------------------------------------------
# Allows another project to write: find_package(Gecko::GeomItf REQUIRED)
install(TARGETS gecko_geom_itf_lib
        EXPORT GeckoGeomItfTargets
        FILE_SET HEADERS DESTINATION include
)

install(EXPORT GeckoGeomItfTargets
        FILE GeckoGeomItfTargets.cmake
        NAMESPACE Gecko::
        DESTINATION lib/cmake/GeckoGeomItf
)
