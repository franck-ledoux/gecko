# -----------------------------------------------------------------------------
# Header-only library
# -----------------------------------------------------------------------------
add_library(gecko_math_lib INTERFACE)
# Alias that allows to write "target_link_libraries(my_exe PRIVATE Gecko::Math)"
add_library(Gecko::Math ALIAS gecko_math_lib)
# Directory where to find header files
target_include_directories(gecko_math_lib INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:include>
)

# List of all the header files to consider
target_sources(gecko_math_lib INTERFACE
        FILE_SET HEADERS
        BASE_DIRS inc
        FILES
        inc/gecko/math/BezierCurve.h
        inc/gecko/math/BezierSurface.h
        inc/gecko/math/BezierVolume.h
        inc/gecko/math/ClosestPoint.h
        inc/gecko/math/CoonsPatch.h
        inc/gecko/math/CurveSurfaceTraits.h
        inc/gecko/math/Point3d.h
        inc/gecko/math/Vector3d.h
)

# Dependency to other Gecko modules
target_link_libraries(gecko_math_lib INTERFACE Gecko::Utils)
# Require C++23
target_compile_features(gecko_math_lib INTERFACE cxx_std_23)

# LIEN IMPORTANT : Injection des flags de couverture sur la bibliothèque
if(TARGET coverage_config)
    target_link_libraries(gecko_math_lib PUBLIC coverage_config)
endif()

# -----------------------------------------------------------------------------
# Test directory
# -----------------------------------------------------------------------------
if (GECKO_BUILD_TESTS)
    add_subdirectory(tst)
endif ()
# -----------------------------------------------------------------------------
# Exportation for external projects
# -----------------------------------------------------------------------------
# Allows another project to write: find_package(Gecko::Math REQUIRED)
install(TARGETS gecko_math_lib
        EXPORT GeckoMathTargets
        FILE_SET HEADERS DESTINATION include
)

install(EXPORT GeckoMathTargets
        FILE GeckoMathTargets.cmake
        NAMESPACE Gecko::
        DESTINATION lib/cmake/GeckoMath
)
