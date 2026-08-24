# -----------------------------------------------------------------------------
# The application layer: the model both front ends drive.
#
# BlockingFacade/GeomModelFacade wrap the kernel in the vocabulary an application
# speaks — ids rather than attribute handles, flat arrays of primitives rather
# than CGAL types, argument checking that throws instead of asserting. Neither is
# specific to Python: the `gecko` extension module and `biy` are both clients, and
# they used to get it by each compiling the same .cpp files out of `python/src`.
# A library, so there is one build of it and the name says what it is.
# -----------------------------------------------------------------------------
add_library(gecko_app_lib STATIC
        src/GeomModelFacade.cpp
        src/BlockingFacade.cpp
)
# Alias that allows to write "target_link_libraries(my_exe PRIVATE Gecko::App)"
add_library(Gecko::App ALIAS gecko_app_lib)

target_include_directories(gecko_app_lib PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:include>
)

target_sources(gecko_app_lib PUBLIC
        FILE_SET HEADERS
        BASE_DIRS inc
        FILES
        inc/gecko/app/BlockingFacade.h
        inc/gecko/app/GeomModelFacade.h
)

target_link_libraries(gecko_app_lib PUBLIC
        Gecko::Block
        Gecko::Geom
        Gecko::Io
        Gecko::Math
        Gecko::Mesh
        Gecko::Utils
)
target_compile_features(gecko_app_lib PUBLIC cxx_std_23)
# Position-independent: the `gecko` Python extension is a shared object and links this in.
set_target_properties(gecko_app_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
