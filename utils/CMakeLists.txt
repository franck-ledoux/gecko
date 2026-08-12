# -----------------------------------------------------------------------------
# Header-only library
# -----------------------------------------------------------------------------
add_library(gecko_utils_lib INTERFACE)
# Alias that allows to write "target_link_libraries(my_exe PRIVATE Gecko::Utils)"
add_library(Gecko::Utils ALIAS gecko_utils_lib)
# Directory where to find header files
target_include_directories(gecko_utils_lib INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:include>
)


# List of all the header files to consider
target_sources(gecko_utils_lib INTERFACE
        FILE_SET HEADERS
        BASE_DIRS inc
        FILES
        inc/gecko/utils/Types.h
        inc/gecko/utils/Groups.h
        inc/gecko/utils/Utils.h
)
# Require C++23
target_compile_features(gecko_utils_lib INTERFACE cxx_std_23)
# -----------------------------------------------------------------------------
# Test directory
# -----------------------------------------------------------------------------
if (GECKO_BUILD_TESTS)
    add_subdirectory(tst)
endif ()
# -----------------------------------------------------------------------------
# Exportation for external projects
# -----------------------------------------------------------------------------
# Allows another project to write: find_package(Gecko::Utils REQUIRED)
install(TARGETS gecko_utils_lib
        EXPORT GeckoUtilsTargets
        FILE_SET HEADERS DESTINATION include
)

install(EXPORT GeckoUtilsTargets
        FILE GeckoUtilsTargets.cmake
        NAMESPACE Gecko::
        DESTINATION lib/cmake/GeckoUtils
)
