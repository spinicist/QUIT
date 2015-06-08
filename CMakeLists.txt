cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(QUIT CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
aux_source_directory(. SRC_LIST)
include(CheckCXXCompilerFlag)

find_package(ITK 4.7 REQUIRED)
include(${ITK_USE_FILE})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(Source)

add_library(quit_core Source/Types.h Source/Util.cpp Source/SignalEquations.cpp Source/Model.cpp Source/Sequence.cpp)
add_custom_target(quit_template SOURCES Source/Filters/ImageToVectorFilter.h Source/Filters/ImageToVectorFilter.hxx
                                        Source/Filters/VectorToImageFilter.h Source/Filters/VectorToImageFilter.hxx
                                        Source/Filters/ApplyAlgorithmFilter.h Source/Filters/ApplyAlgorithmFilter.hxx
                                        Source/Filters/ReorderVectorFilter.h Source/Filters/ReorderVectorFilter.hxx
                                        Source/RegionContraction.h)

set(PROGRAMS qnewimage qcomplex qsignal
             qdespot1 qdespot2 qdespot1hifi qdespot2fm qmcdespot
             qafi qssfpbands
             qmultiecho)

foreach(PROGRAM ${PROGRAMS})
	add_executable(${PROGRAM} Source/${PROGRAM}.cpp)
	target_link_libraries(${PROGRAM} quit_core ${ITK_LIBRARIES})
endforeach(PROGRAM)

set(TESTS test_common.sh test_scdespot.sh test_mcdespot.sh test_relax.sh test_all.sh)
foreach(TEST ${TESTS})
	configure_file(Test/${TEST} Test/${TEST} COPYONLY)
endforeach(TEST)
