project(QUITK)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)

set(CMAKE_BUILD_TYPE Release)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
	message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

find_package(ITK 4.7 REQUIRED)
include(${ITK_USE_FILE})

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(despot1 Source/despot1.cpp Source/SignalEquations.cpp Source/Model.cpp Source/Algorithm.cpp Source/Sequence.cpp Source/Util.cpp Source/Filters/ImageToVectorFilter.h Source/Filters/ImageToVectorFilter.hxx Source/Filters/DESPOT1Filter.h Source/Filters/DESPOT1Filter.hxx)
target_link_libraries(despot1 ${ITK_LIBRARIES})
