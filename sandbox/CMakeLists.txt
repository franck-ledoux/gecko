# -----------------------------------------------------------------------------
# Scratch executable for local development. Not installed, not part of the public API/docs (see
# GECKO_BUILD_SANDBOX in the root CMakeLists.txt) — edit sandbox/src/main.cpp freely.
# -----------------------------------------------------------------------------
add_executable(gecko_sandbox src/main.cpp)
target_link_libraries(gecko_sandbox PRIVATE
        Gecko::Block
        Gecko::Geom
        Gecko::Io
)
target_compile_features(gecko_sandbox PRIVATE cxx_std_23)
