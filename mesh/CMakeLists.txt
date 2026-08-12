# -----------------------------------------------------------------------------
# Header-only library
# -----------------------------------------------------------------------------
add_library(gecko_mesh_lib INTERFACE)
# Alias that allows to write "target_link_libraries(mon_exe PRIVATE Gecko::Math)"
add_library(Gecko::Mesh ALIAS gecko_mesh_lib)

# Directory where to find header files
target_include_directories(gecko_mesh_lib INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:include>
)
# List of all the header files to consider
target_sources(gecko_mesh_lib INTERFACE
        FILE_SET HEADERS
        BASE_DIRS inc
        FILES
        inc/gecko/mesh/UnstructuredMesh.h
        inc/gecko/mesh/VariableRegistry.h
)
# Dependency to other Gecko modules
target_link_libraries(gecko_mesh_lib INTERFACE Gecko::Math Gecko::Utils)
# Require C++23
target_compile_features(gecko_mesh_lib INTERFACE cxx_std_23)
# -----------------------------------------------------------------------------
# Test directory
# -----------------------------------------------------------------------------
if (GECKO_BUILD_TESTS)
    add_subdirectory(tst)
endif ()
# -----------------------------------------------------------------------------
# Exportation for external projects
# -----------------------------------------------------------------------------
# Allows another project to write: find_package(GeckoMesh REQUIRED)
install(TARGETS gecko_mesh_lib
        EXPORT GeckoMeshTargets
        FILE_SET HEADERS DESTINATION include
)

install(EXPORT GeckoMeshTargets
        FILE GeckoMeshTargets.cmake
        NAMESPACE Gecko::
        DESTINATION lib/cmake/GeckoMesh
)
