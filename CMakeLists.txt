cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(QUIT CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_FLAGS -fpermissive)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
aux_source_directory(. SRC_LIST)
include(CheckCXXCompilerFlag)
include(CMakeToolsHelpers OPTIONAL)

# Prevent undefined lookup errors when building libqi as a shared library with Clang
message( STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}" )
if(CMAKE_CXX_COMPILER_ID MATCHES "(C|c?)lang")
    set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

find_package(ITK 4.10.0 REQUIRED)
include(${ITK_USE_FILE})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories( ${PROJECT_SOURCE_DIR}/Source )
include_directories( ${PROJECT_SOURCE_DIR}/External/CppNumericalSolvers/include )

# Setup version
include_directories(${CMAKE_CURRENT_BINARY_DIR}) # For version file
set( VERSION_FILE_NAME "VersionFile")
add_custom_target( qi_version ALL DEPENDS version_depend_dummy )
add_custom_command(OUTPUT version_depend_dummy
                   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/CMake/GetVersion.cmake)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${VERSION_FILE_NAME} PROPERTIES GENERATED TRUE HEADER_FILE_ONLY TRUE )

set(QIDIR "Source/QI")
set(QI "qi")
add_library( ${QI} STATIC
             ${QIDIR}/Types.h ${QIDIR}/Util.cpp ${QIDIR}/ThreadPool.cpp
             ${QIDIR}/Option.cpp ${QIDIR}/Args.cpp
             ${QIDIR}/Kernels.cpp ${QIDIR}/Lineshape.cpp
             ${QIDIR}/GoldenSection.cpp )
add_dependencies( ${QI} qi_version )
add_custom_target(qi_template SOURCES Source/Filters/ImageToVectorFilter.h Source/Filters/ImageToVectorFilter.hxx
                                      Source/Filters/VectorToImageFilter.h Source/Filters/VectorToImageFilter.hxx
                                      Source/Filters/ApplyAlgorithmFilter.h Source/Filters/ApplyAlgorithmFilter.hxx
                                      Source/Filters/ReorderVectorFilter.h Source/Filters/ReorderVectorFilter.hxx
                                      Source/Filters/ReorderImageFilter.h Source/Filters/ReorderImageFilter.hxx
                                      ${QIDIR}/RegionContraction.h)

set( SIGDIR ${QIDIR}/Signals )
set( SIGLIB "qisig" )
add_library( ${SIGLIB} STATIC
             ${SIGDIR}/Common.cpp ${SIGDIR}/SignalEquations.cpp
             ${SIGDIR}/SPGR.cpp ${SIGDIR}/SSFP.cpp ${SIGDIR}/SSFP_MC.cpp ${SIGDIR}/MPRAGE.cpp )

set( SEQDIR ${QIDIR}/Sequences )
set( SEQLIB "qiseq" )
add_library( ${SEQLIB} STATIC
             ${SEQDIR}/SequenceBase.cpp ${SEQDIR}/Sequences.cpp
             ${SEQDIR}/SteadyStateSequence.cpp ${SEQDIR}/MPRAGESequence.cpp ${SEQDIR}/SpinEcho.cpp )

set( MDLDIR ${QIDIR}/Models )
set( MDLLIB "qimdl" )
add_library( ${MDLLIB} STATIC
             ${MDLDIR}/Model.cpp ${MDLDIR}/DESPOT_2C.cpp ${MDLDIR}/DESPOT_3C.cpp ${MDLDIR}/qMT.cpp )

set( ALGDIR ${QIDIR}/Algorithms )
set( ALGLIB "qialg" )
add_library( ${ALGLIB} STATIC
             ${ALGDIR}/Banding.cpp ${ALGDIR}/Ellipse.cpp ${ALGDIR}/EllipseHelpers.cpp ${ALGDIR}/EllipseFit.cpp
             ${ALGDIR}/CEST.cpp )

set( UTILS qihdr qinewimage qidiff 
           qicomplex qiaffine qimask qikfilter qireorder
           qisplitsubjects qiunwrap qipolyfit qipolyimg
           qi_roivolumes qi_roivalues
           qi_glmsetup qi_glmcontrasts )
set( QUANT qisignal qidespot1 qidespot2 qidespot1hifi qidespot2fm qimcdespot qiafi qimultiecho qimp2rage qicest )
set( ELLIPSE qissfpbands qiesmap )

foreach(PROGRAM ${UTILS})
    add_executable(${PROGRAM} Source/${PROGRAM}.cpp)
    target_link_libraries(${PROGRAM} ${QI} ${ITK_LIBRARIES})
endforeach(PROGRAM)

foreach(PROGRAM ${QUANT})
    add_executable(${PROGRAM} Source/${PROGRAM}.cpp)
    target_link_libraries(${PROGRAM} ${ALGLIB} ${MDLLIB} ${SEQLIB} ${SIGLIB} ${QI} ${ITK_LIBRARIES})
endforeach(PROGRAM)

foreach(PROGRAM ${ELLIPSE})
    add_executable(${PROGRAM} Source/${PROGRAM}.cpp)
    target_link_libraries(${PROGRAM} ${ALGLIB} ${MDLLIB} ${SEQLIB} ${SIGLIB} ${QI} ${ITK_LIBRARIES})
endforeach(PROGRAM)

set(TESTS_DIR Test)
set(TESTS despot1.bats test_common.sh test_d1.sh test_fm.sh test_2C.sh test_3C.sh test_es.sh
          test_relax.sh test_ssfp.sh test_poly.sh test_all.sh)
foreach(TEST ${TESTS})
	configure_file(${TESTS_DIR}/${TEST} ${TESTS_DIR}/${TEST} COPYONLY)
endforeach(TEST)

enable_testing()
add_test( NAME DESPOT1 WORKING_DIRECTORY ${TESTS_DIR} 
          COMMAND ${PROJECT_SOURCE_DIR}/External/bats/bin/bats despot1.bats )

set(SCRIPTS_DIR Scripts)
set(SCRIPTS example_mcd_hifi.sh example_mcd_b1.sh)
foreach(SCRIPT ${SCRIPTS})
    configure_file(${SCRIPTS_DIR}/${SCRIPT} ${SCRIPTS_DIR}/${SCRIPT} COPYONLY)
endforeach(SCRIPT)

set(PROGRAMS ${UTILS} ${QUANT} ${ELLIPSE} )
install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )
