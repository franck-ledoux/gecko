# -----------------------------------------------------------------------------
# Header-only library
# -----------------------------------------------------------------------------
add_library(gecko_io_lib INTERFACE)
# Alias that allows to write "target_link_libraries(my_exe PRIVATE Gecko::Io)"
add_library(Gecko::Io ALIAS gecko_io_lib)

# Directory where to find header files
target_include_directories(gecko_io_lib INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:include>
)
# List of all the header files to consider
target_sources(gecko_io_lib INTERFACE
        FILE_SET HEADERS
        BASE_DIRS inc
        FILES
        inc/gecko/io/GmshCommon.h
        inc/gecko/io/GmshMeshReader.h
        inc/gecko/io/GmshMeshWriter.h
        inc/gecko/io/VtkMeshWriter.h
)
# Dependency to other Gecko modules
target_link_libraries(gecko_io_lib INTERFACE Gecko::Mesh Gecko::Utils)
# Require C++23
target_compile_features(gecko_io_lib INTERFACE cxx_std_23)
# -----------------------------------------------------------------------------
# Test directory
# -----------------------------------------------------------------------------
if (GECKO_BUILD_TESTS)
    add_subdirectory(tst)
endif ()
# -----------------------------------------------------------------------------
# Exportation for external projects
# -----------------------------------------------------------------------------
# Allows another project to write: find_package(Gecko::Io REQUIRED)
install(TARGETS gecko_io_lib
        EXPORT GeckoIoTargets
        FILE_SET HEADERS DESTINATION include
)

install(EXPORT GeckoIoTargets
        FILE GeckoIoTargets.cmake
        NAMESPACE Gecko::
        DESTINATION lib/cmake/GeckoIo
)
