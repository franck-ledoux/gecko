# Install script for directory: /Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Library/Developer/CommandLineTools/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/AdolcForward"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/AlignedVector3"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/ArpackSupport"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/AutoDiff"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/BVH"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/EulerAngles"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/FFT"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/IterativeSolvers"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/KroneckerProduct"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/LevenbergMarquardt"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/MatrixFunctions"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/MoreVectorization"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/MPRealSupport"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/NonLinearOptimization"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/NumericalDiff"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/OpenGLSupport"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/Polynomials"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/Skyline"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/SparseExtra"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/SpecialFunctions"
    "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/Splines"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/Users/bibi/Developpement/gecko/build/_deps/eigen3-src/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/bibi/Developpement/gecko/build/_deps/eigen3-build/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/bibi/Developpement/gecko/build/_deps/eigen3-build/unsupported/Eigen/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
