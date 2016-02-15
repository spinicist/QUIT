cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(QUIT CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_FLAGS -fpermissive)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
aux_source_directory(. SRC_LIST)
include(CheckCXXCompilerFlag)

find_package(ITK 4.8.1 REQUIRED)
include(${ITK_USE_FILE})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(Source)

# Setup version
include_directories(${CMAKE_CURRENT_BINARY_DIR}) # For version file
add_custom_target( qi_version ALL
                              DEPENDS VERSION_DEPENDENCY )
message( STATUS "Source: ${CMAKE_CURRENT_SOURCE_DIR}" )
add_custom_command(OUTPUT VERSION_DEPENDENCY
                   COMMAND ${CMAKE_COMMAND} -DCMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/CMake/GetVersion.cmake)

set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/BuildVersion.h PROPERTIES GENERATED TRUE HEADER_FILE_ONLY TRUE )

add_library(qi_core Source/Types.h Source/Util.cpp Source/ThreadPool.cpp Source/SignalEquations.cpp Source/Model.cpp Source/Sequence.cpp)
add_dependencies( qi_core qi_version )
add_custom_target(qi_template SOURCES Source/Filters/ImageToVectorFilter.h Source/Filters/ImageToVectorFilter.hxx
                                      Source/Filters/VectorToImageFilter.h Source/Filters/VectorToImageFilter.hxx
                                      Source/Filters/ApplyAlgorithmFilter.h Source/Filters/ApplyAlgorithmFilter.hxx
                                      Source/Filters/ReorderVectorFilter.h Source/Filters/ReorderVectorFilter.hxx
                                      Source/Filters/ReorderImageFilter.h Source/Filters/ReorderImageFilter.hxx
                                      Source/RegionContraction.h)

set(PROGRAMS qihdr qinewimage qicomplex qisignal qiaffine qimask qikfilter qireorder qisplitsubjects
             qidespot1 qidespot2 qidespot1hifi qidespot2fm qimcdespot
             qiafi qissfpbands qiunwrap
             qimultiecho qimp2rage )

foreach(PROGRAM ${PROGRAMS})
	add_executable(${PROGRAM} Source/${PROGRAM}.cpp)
	target_link_libraries(${PROGRAM} qi_core ${ITK_LIBRARIES})
endforeach(PROGRAM)

set(TESTS_DIR Test)
set(TESTS test_common.sh test_d1.sh test_fm.sh test_2C.sh test_3C.sh test_relax.sh test_ssfp.sh test_all.sh)
foreach(TEST ${TESTS})
	configure_file(${TESTS_DIR}/${TEST} ${TESTS_DIR}/${TEST} COPYONLY)
endforeach(TEST)

set(SCRIPTS_DIR Scripts)
set(SCRIPTS example_mcd_hifi.sh example_mcd_b1.sh example_mcd_fsl_registration.sh)
foreach(SCRIPT ${SCRIPTS})
	configure_file(${SCRIPTS_DIR}/${SCRIPT} ${SCRIPTS_DIR}/${SCRIPT} COPYONLY)
endforeach(SCRIPT)

install(TARGETS ${PROGRAMS} RUNTIME DESTINATION bin)
