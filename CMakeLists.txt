project(QUITK)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)

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

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")

find_package(ITK 4.7 REQUIRED)
include(${ITK_USE_FILE})

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_library(quitk_core Source/SignalEquations.cpp Source/Model.cpp Source/Sequence.cpp Source/Util.cpp)
add_custom_target(quitk_templ SOURCES Source/Filters/ImageToVectorFilter.h Source/Filters/ImageToVectorFilter.hxx
                                      Source/Filters/VectorToImageFilter.h Source/Filters/VectorToImageFilter.hxx
                                      Source/Filters/ApplyAlgorithmFilter.h Source/Filters/ApplyAlgorithmFilter.hxx)

set(PROGRAMS despot1 despot1hifi despot2 mcsignal ssfpgs)
foreach(PROGRAM ${PROGRAMS})
	add_executable(${PROGRAM} Source/${PROGRAM}.cpp)
	target_link_libraries(${PROGRAM} quitk_core ${ITK_LIBRARIES})
endforeach(PROGRAM)

set(TESTS test_common.sh test_scdespot.sh test_all.sh)
foreach(TEST ${TESTS})
	configure_file(Test/${TEST} Test/${TEST} COPYONLY)
endforeach(TEST)
