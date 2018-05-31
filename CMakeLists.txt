cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(QUIT CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)
aux_source_directory(. SRC_LIST)
include(CheckCXXCompilerFlag)
include(CMakeToolsHelpers OPTIONAL)

set( CMAKE_CXX_FLAGS "-Wall -Wextra" )
set( Args_DIR "${PROJECT_SOURCE_DIR}/External/Args"
     CACHE PATH "Path to Args library (usually External/Args" )
set( Cereal_DIR "${PROJECT_SOURCE_DIR}/External/cereal/include"
     CACHE PATH "Path to Cereal library (usually External/Cereal/include" )
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)
include( ${PROJECT_SOURCE_DIR}/CMake/BuildType.cmake )

if( NOT DEFINED EIGEN3_INCLUDE_DIRS )
    find_package( Eigen3 REQUIRED )
endif()
message( STATUS "Found Eigen in ${EIGEN3_INCLUDE_DIRS}" )
include_directories( SYSTEM ${EIGEN3_INCLUDE_DIRS} )

if( NOT DEFINED CERES_INCLUDE_DIRS )
    find_package( Ceres 1.12.0 )
endif()
include_directories( SYSTEM ${CERES_INCLUDE_DIRS} )

find_package(ITK 4.10.0 REQUIRED
              COMPONENTS
                ITKBinaryMathematicalMorphology ITKCommon ITKConnectedComponents
                ITKFFT ITKIOImageBase ITKIONRRD ITKIONIFTI
                ITKIOTransformBase ITKImageCompose
                ITKImageFeature ITKImageFilterBase ITKImageFunction ITKImageGrid
                ITKImageIntensity ITKImageStatistics ITKLabelMap ITKLabelVoting
                ITKMathematicalMorphology ITKOptimizers ITKRegistrationCommon
                ITKSmoothing ITKThresholding ITKTransform )
include( ${ITK_USE_FILE} )
include_directories( SYSTEM ${Args_DIR} ${Cereal_DIR} )

OPTION(BUILD_SHARED_LIBS "Build shared libraries." OFF)

set( VERSION_MAJOR "2" )
set( VERSION_MINOR "0" )
set( VERSION_PATCH "1" )

add_subdirectory( Source )
add_subdirectory( Scripts )

include( InstallRequiredSystemLibraries )
set( CPACK_PACKAGE_DESCRIPTION_SUMMARY "Quantitative Imaging Tools" )
set( CPACK_PACKAGE_CONTACT "tobias.wood@kcl.ac.uk" )
set( CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md" )
set( CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR} )
set( CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR} )
set( CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH} )
set( CPACK_PACKAGING_INSTALL_PREFIX "" )
set( CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt" )
if( APPLE )
    set( CPACK_GENERATOR TGZ )
elseif( UNIX AND NOT APPLE )
    set( CPACK_GENERATOR TGZ DEB )
else()
    message( STATUS "Using CPack default generators (all of them)" )
endif()
include( CPack )
