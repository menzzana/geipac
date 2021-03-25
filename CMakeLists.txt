project(geipac)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

###############################################################################
# Link to static libraries
###############################################################################

if (STATIC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
  set(Boost_USE_STATIC_LIBS ON)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif(STATIC)

###############################################################################
# Check for dependencies
###############################################################################

find_package(Boost COMPONENTS filesystem system program_options REQUIRED)
if(Boost_FOUND)
  message(STATUS "Package Boost-${Boost_VERSION} found : ${Boost_INCLUDE_DIR}")
else(Boost_FOUND)
  message(FATAL_ERROR "The package Boost has not been found")
endif(Boost_FOUND)

set(EIGEN3_INCLUDE_DIR ${EIGEN_INCLUDE_DIR})
if(EIGEN3_INCLUDE_DIR)
  message(STATUS "Eigen3 found")
else(EIGEN3_INCLUDE_DIR)
  message(FATAL_ERROR "Eigen3 installation has not been found")
endif(EIGEN3_INCLUDE_DIR)

###############################################################################
# if OpenMP not found the code will be compiled for serial execution
# Comment the find of OpenMP to compile and run serially
###############################################################################

find_package(OpenMP)
if (OPENMP_FOUND)
  message(STATUS "Package OpenMP found.")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else(OPENMP_FOUND)
  message(STATUS "The package OpenMP was not found.")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSERIAL")
endif(OPENMP_FOUND)

###############################################################################
# Link and include
###############################################################################
include(version)
add_executable(${PROJECT_NAME} ${SRC_LIST})
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
link_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
