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
set( Cereal_DIR "${PROJECT_SOURCE_DIR}/External/cereal/include"
     CACHE PATH "Path to Cereal library (usually External/Cereal/include" )
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)
include( ${PROJECT_SOURCE_DIR}/CMake/BuildType.cmake )
find_package(Eigen3 REQUIRED)
find_package(Ceres 1.12.0 REQUIRED COMPONENTS C++11)
find_package(ITK 4.10.0 REQUIRED)

include_directories( ${EIGEN3_INCLUDE_DIR} )
include_directories( ${CERES_INCLUDE_DIRS} )
include( ${ITK_USE_FILE} )
include_directories( ${Args_DIR} ${Cereal_DIR} )

OPTION(BUILD_SHARED_LIBS "Build shared libraries." ON)

add_subdirectory( Source )
add_subdirectory( Scripts )

INCLUDE(InstallRequiredSystemLibraries)
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Quantitative Imaging Tools")
SET(CPACK_PACKAGE_VENDOR "tobias.wood@kcl.ac.uk")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
SET(CPACK_PACKAGE_VERSION_MAJOR "2")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_GENERATOR "STGZ")
set(CPACK_SET_DESTDIR true)
INCLUDE(CPack)