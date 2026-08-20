# -----------------------------------------------------------------------------
# biy ("Block It Yourself"): interactive 3D viewer for building a Blocking against a GeomModel.
# Polyscope drives the render loop on the main thread (via its non-blocking frameTick()), while an
# embedded Python interpreter runs a REPL on a background thread — both acting on the very same
# facade objects. Behind GECKO_BUILD_BIY (root CMakeLists.txt); see docs/user-guide/biy.md.
# -----------------------------------------------------------------------------
find_package(Python COMPONENTS Interpreter Development.Embed REQUIRED)

add_executable(biy
        src/main.cpp
        src/BiyApp.cpp
        src/BiyConfig.cpp
        src/PythonConsole.cpp
        # The façade + binding glue, shared verbatim with the `gecko` Python extension module. Too
        # small to be worth an installed library of its own (same call as python/CMakeLists.txt).
        ${CMAKE_SOURCE_DIR}/python/src/GeomModelFacade.cpp
        ${CMAKE_SOURCE_DIR}/python/src/BlockingFacade.cpp
        ${CMAKE_SOURCE_DIR}/python/src/PyGeckoBindings.cpp
)
target_include_directories(biy PRIVATE ${CMAKE_SOURCE_DIR}/python/src)
target_link_libraries(biy PRIVATE
        Gecko::Block
        Gecko::Geom
        Gecko::Io
        polyscope
        pybind11::embed
        json_dep
)
target_compile_features(biy PRIVATE cxx_std_23)
