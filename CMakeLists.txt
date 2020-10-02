cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
cmake_policy(VERSION 3.12)
project(QUIT CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/CMake)
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${PROJECT_SOURCE_DIR}/CMake)
include(CheckCXXCompilerFlag)
include(CMakeToolsHelpers OPTIONAL)
include(CMakePrintHelpers)

include( ${PROJECT_SOURCE_DIR}/CMake/BuildType.cmake )
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++" CACHE INTERNAL "" FORCE)
endif()

# Configure directory paths
set( External_Include_DIR "${PROJECT_SOURCE_DIR}/External/include"
     CACHE PATH "Path to external library includes (usually External/include" )
include_directories( SYSTEM ${External_Include_DIR} )
find_package(Eigen3 CONFIG REQUIRED)
find_package(Ceres CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(args CONFIG REQUIRED)
find_package(ITK 5.0.0 REQUIRED
              COMPONENTS
                ITKBinaryMathematicalMorphology
                ITKCommon
                ITKConnectedComponents
                ITKFFT
                ITKImageCompose
                ITKImageFeature
                ITKImageFilterBase
                ITKImageIntensity
                ITKImageStatistics
                ITKLabelMap
                ITKLabelVoting
                ITKMathematicalMorphology
                ITKThresholding
                ITKIOTransformInsightLegacy
                ITKIONIFTI )
include( ${ITK_USE_FILE} )

add_subdirectory( Source )
