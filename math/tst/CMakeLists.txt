add_executable(test_math
        bezier_curve_tests.cpp
        bezier_surface_tests.cpp
        bezier_volume_tests.cpp
        closest_point_tests.cpp
        coons_patch_tests.cpp
        PointVector_tests.cpp)
target_link_libraries(test_math PRIVATE
    Gecko::Math
    Catch2::Catch2
    Catch2::Catch2WithMain
)
# To get access to the generated file "unit_test_config.h"
target_include_directories(test_math PRIVATE ${CMAKE_BINARY_DIR})

# LIEN IMPORTANT : Injection des flags de couverture sur l'exécutable de tests
if(TARGET coverage_config)
    target_link_libraries(test_math PRIVATE coverage_config)
endif()
add_test(NAME test_math COMMAND test_math)
