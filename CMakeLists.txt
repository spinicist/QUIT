cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
cmake_policy(VERSION 3.12)
project(QUIT CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)

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
if(${CMAKE_SYSTEM_VERSION} EQUAL "Linux")
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
endif()

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

add_subdirectory( Source )
