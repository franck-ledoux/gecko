add_library(mcts_c_lib
    src/BlockingAction.cpp
    src/BlockingState.cpp
    src/BlockingRewardFunction.cpp
    spam/src/MCTSAgent.cpp
    spam/src/MCTSLoop.cpp
    spam/src/MCTSSelectionFunction.cpp
    spam/src/MCTSTree.cpp
)
target_include_directories(mcts_c_lib PUBLIC inc spam/inc)
target_link_libraries(mcts_c_lib PUBLIC
    cgal_dep
    boost_dep
    nlohmann_json::nlohmann_json
    gmds_core_cad_lib
    gmds_core_cadfac_lib
    gmds_core_ig_lib
    gmds_core_igalgo_lib
    gmds_core_io_lib
    gmds_core_math_lib
    gmds_core_utils_lib
    cblock_lib
    Eigen3::Eigen
)

add_executable(mcts_c_blocker src/main_blocker.cpp)
target_link_libraries(mcts_c_blocker PRIVATE mcts_c_lib)
