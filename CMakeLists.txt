project(geipac)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11 -fno-delete-null-pointer-checks")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

###############################################################################
# Build static libraries
###############################################################################
if (BUILD_STATIC_EXE)
  message(STATUS "Geipac will be built with static libraries.")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
  set(Boost_USE_STATIC_LIBS ON)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)
  set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
  set(CMAKE_SHARED_LIBRARY_C_FLAGS)
  set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
  set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
  set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
else (BUILD_STATIC_EXE)
  message(STATUS "Geipac will be built with dynamic libraries.")
endif(BUILD_STATIC_EXE)

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
# if OpenMP not found or serial was selected 
# the code will be compiled for serial execution
###############################################################################

find_package(OpenMP)
if (NOT OPENMP_FOUND OR BUILD_SERIAL_EXE)
  message(STATUS "Geipac will be built for serial execution.")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSERIAL")  
else (NOT OPENMP_FOUND OR BUILD_SERIAL_EXE)
  message(STATUS "Geipac will be built with OpenMP for parallel execution.")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(NOT OPENMP_FOUND OR BUILD_SERIAL_EXE)

###############################################################################
# Set version
###############################################################################

set(VERSION_MAJOR "0")
set(VERSION_MINOR "3")
set(VERSION_PATCH "1")
include(version)
message(STATUS "Building Geipac V${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}-${GIT_HASH}")

###############################################################################
# Link and include
###############################################################################

add_executable(${PROJECT_NAME} ${SRC_LIST})
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
link_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
