cmake_minimum_required(VERSION 3.10.2 FATAL_ERROR)

project(QUIT CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)
add_compile_options("-Werror" "-Wall" "-Wextra" "-Wshadow" "-ftemplate-backtrace-limit=30")
set(ASAN $<$<CONFIG:RelWithDebInfo>:-fsanitize=address>)
add_compile_options(${ASAN})
add_link_options(${ASAN})
aux_source_directory(. SRC_LIST)
include(CheckCXXCompilerFlag)
include(CMakeToolsHelpers OPTIONAL)

include( ${PROJECT_SOURCE_DIR}/CMake/BuildType.cmake )
set( External_Include_DIR "${PROJECT_SOURCE_DIR}/External/include"
     CACHE PATH "Path to external library includes (usually External/include" )
include_directories( SYSTEM ${External_Include_DIR} )

# Configure directory paths
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/CMake)
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${PROJECT_SOURCE_DIR}/CMake)

if( NOT DEFINED EIGEN3_INCLUDE_DIR )
    find_package( Eigen3 REQUIRED )
endif()
message( STATUS "Found Eigen in ${EIGEN3_INCLUDE_DIR}" )
include_directories( SYSTEM ${EIGEN3_INCLUDE_DIR} )

if(DEFINED CERES_BUILD_DIR)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CERES_BUILD_DIR})
endif()
find_package( Ceres 1.14.0 )
if(CERES_FOUND)
    include_directories( SYSTEM ${CERES_INCLUDE_DIRS} )
endif()

# fmtlib: Sane output
find_package(fmtlib REQUIRED CONFIG)
if(fmtlib_FOUND)
    add_library(fmtlib STATIC ${FMTLIB_SOURCES})
    target_include_directories(fmtlib SYSTEM PUBLIC ${FMTLIB_INCLUDE_DIR})
endif()

find_package(ITK 5.0.0 REQUIRED
              COMPONENTS
                ITKBinaryMathematicalMorphology
                ITKCommon
                ITKConnectedComponents
                ITKFFT
                ITKIOImageBase
                ITKIOTransformBase
                ITKImageCompose
                ITKImageFeature
                ITKImageFilterBase
                ITKImageGrid
                ITKImageIntensity
                ITKImageStatistics
                ITKLabelMap
                ITKLabelVoting
                ITKMathematicalMorphology
                ITKThresholding
                ITKTransform
                ITKImageIO
                ITKTransformIO
                ITKIONIFTI )
include( ${ITK_USE_FILE} )

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
