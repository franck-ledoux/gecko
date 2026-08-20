add_executable(test_block
        blocking_classification_tests.cpp
        blocking_connectivity_tests.cpp
        blocking_creation_tests.cpp
        blocking_curve_evaluation_tests.cpp
        blocking_edit_tests.cpp
        blocking_mesh_generation_tests.cpp
        blocking_node_move_tests.cpp
        cell_data_tests.cpp)
target_link_libraries(test_block PRIVATE
    Gecko::Block
    Gecko::Geom
    Catch2::Catch2
    Catch2::Catch2WithMain
)
# To get access to the generated file "unit_test_config.h"
target_include_directories(test_block PRIVATE ${CMAKE_BINARY_DIR})

# LIEN IMPORTANT : Injection des flags de couverture sur l'exécutable de tests
if(TARGET coverage_config)
    target_link_libraries(test_block PRIVATE coverage_config)
endif()

add_test(NAME test_block COMMAND test_block)
