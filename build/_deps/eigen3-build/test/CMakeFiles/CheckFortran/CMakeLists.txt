cmake_minimum_required(VERSION 3.31.5)
set(CMAKE_MODULE_PATH "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/cmake")
project(CheckFortran Fortran)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/result.cmake"
  "set(CMAKE_Fortran_COMPILER \"${CMAKE_Fortran_COMPILER}\")\n"
  
  )