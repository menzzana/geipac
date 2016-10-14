project(geipac)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11")

###############################################################################
# Set version
###############################################################################

set(GEIPAC_VERSION_MAJOR "0")
set(GEIPAC_VERSION_MINOR "1")
set(GEIPAC_VERSION_PATCH "0")
configure_file( version_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/version_config.h )

###############################################################################
# Check for dependencies
# if mpi not found parma will be compiled for serial execution
# Uncomment the find of MPI to compile and run parma serially
###############################################################################

find_package(Boost COMPONENTS filesystem system program_options REQUIRED)
if(Boost_FOUND)
  message(STATUS "Package Boost-${Boost_VERSION} found : ${Boost_INCLUDE_DIR}")
else(Boost_FOUND)
  message(FATAL_ERROR "The package Boost has not been found")
endif(Boost_FOUND)
link_directories(${Boost_INCLUDE_DIRS})

set(EIGEN3_INCLUDE_DIR ${EIGEN_INCLUDE_DIR})
if(EIGEN3_INCLUDE_DIR)
  message(STATUS "Eigen3 found")
else(EIGEN3_INCLUDE_DIR)
  message(FATAL_ERROR "Eigen3 installation has not been found")
endif(EIGEN3_INCLUDE_DIR)

find_package(MPI)
if(MPI_FOUND)
  message(STATUS "Package MPI found.")
else(MPI_FOUND)
  message(STATUS "The package MPI was not found.")
  set(CMAKE_CXX_FLAGS "-DSERIAL")
endif(MPI_FOUND)

###############################################################################
# Link and include
###############################################################################

include_directories(${CMAKE_CURRENT_BINARY_DIR} ${MPI_CXX_INCLUDE_PATH} ${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
