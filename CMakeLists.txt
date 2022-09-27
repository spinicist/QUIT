cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
cmake_policy(VERSION 3.12)
project(QUIT CXX)

# Configure directory paths
set( External_Include_DIR "${PROJECT_SOURCE_DIR}/External/include"
     CACHE PATH "Path to external library includes (usually External/include" )
include_directories( SYSTEM ${External_Include_DIR} )
find_package(Eigen3 REQUIRED)
find_package(Ceres REQUIRED)
find_package(fmt REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(args REQUIRED)
find_package(ITK 5.3.0 REQUIRED
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
                ITKIONIFTI
             )
include( ${ITK_USE_FILE} )

add_subdirectory( Source )
