cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(QUIT CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)
aux_source_directory(. SRC_LIST)
include(CheckCXXCompilerFlag)
include(CMakeToolsHelpers OPTIONAL)

# Prevent undefined lookup errors when building libqi as a shared library with Clang
message( STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}" )
if(CMAKE_CXX_COMPILER_ID MATCHES "(C|c?)lang")
    set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

set( Args_DIR "${PROJECT_SOURCE_DIR}/External/Args"
     CACHE PATH "Path to Args library (usually External/Args" )
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)
find_package(Eigen3 REQUIRED)
find_package(Ceres 1.12.0 REQUIRED COMPONENTS C++11)
find_package(ITK 4.10.0 REQUIRED)

include_directories( ${PROJECT_SOURCE_DIR}/Source )
include_directories( ${EIGEN3_INCLUDE_DIR} )
include_directories( ${CERES_INCLUDE_DIRS} )
include( ${ITK_USE_FILE} )
include_directories( ${Args_DIR} )

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
             ${QIDIR}/Option.cpp
             ${QIDIR}/Kernels.cpp ${QIDIR}/Lineshape.cpp ${QIDIR}/Fit.cpp )
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

set( ELDIR ${QIDIR}/Ellipse )
add_library( CEST    STATIC ${QIDIR}/CEST.cpp )
target_link_libraries( CEST ${CERES_LIBRARIES} )
add_library( Banding STATIC ${QIDIR}/Banding.cpp )
add_library( Ellipse STATIC ${ELDIR}/EllipseFilters.cpp ${ELDIR}/EllipseHelpers.cpp ${ELDIR}/EllipseFits.cpp )
target_link_libraries( Ellipse ${CERES_LIBRARIES} ${SIGLIB} )
add_library( Masking STATIC ${QIDIR}/Masking.cpp )

set( UTILS qihdr qinewimage qidiff qiconvert
           qicomplex qiaffine qimask qikfilter qireorder
           qisplitsubjects qiunwrap qipolyfit qipolyimg
           qi_coil_combine qi_bet
           qi_rois qi_glmsetup qi_glmcontrasts )
set( QUANT qisignal
           qiafi qidream
           qidespot1 qidespot2 qidespot1hifi qidespot2fm
           qimcdespot 
           qimultiecho qimp2rage
           qissfpbands qi_ellipse qi_emt qicest )

foreach(PROGRAM ${UTILS})
    add_executable(${PROGRAM} Source/${PROGRAM}.cpp)
    target_link_libraries(${PROGRAM} ${QI} ${ITK_LIBRARIES})
endforeach(PROGRAM)

foreach(PROGRAM ${QUANT})
    add_executable(${PROGRAM} Source/${PROGRAM}.cpp)
    target_link_libraries(${PROGRAM} ${MDLLIB} ${SEQLIB} ${SIGLIB} ${QI} 
                          ${ITK_LIBRARIES} ${CERES_LIBRARIES})
endforeach(PROGRAM)

## Add extra libraries here
target_link_libraries( qicest          CEST )
target_link_libraries( qissfpbands     Banding )
target_link_libraries( qi_ellipse      Ellipse )
target_link_libraries( qisplitsubjects Masking )
target_link_libraries( qimask          Masking )

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
set(SCRIPTS qi_composer.sh)
foreach(SCRIPT ${SCRIPTS})
    INSTALL( FILES ${SCRIPTS_DIR}/${SCRIPT} 
             DESTINATION bin 
             PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
                         OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE )
endforeach(SCRIPT)

set(PROGRAMS ${UTILS} ${QUANT} ${ELLIPSE} )
install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )
