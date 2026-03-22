add_library(gmds_core_utils_lib
    src/BitVector.cpp
    src/CommonTypes.cpp
    src/VariableManager.cpp
)
target_include_directories(gmds_core_utils_lib PUBLIC inc)

add_subdirectory(tst)
