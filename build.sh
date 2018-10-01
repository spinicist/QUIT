#!/bin/bash -eux

#
# easy_build.sh
#
# A hopefully simple script to download all dependencies for QUIT and then
# compile automatically, for those who are unfamiliar with CMake etc.
#
# Requires:
#
# 1. CMake version 3.2 or greater
# 2. A C++11 compliant compiler e.g. GCC 4.8 or higher
#

WD=$PWD

# Command line options
CERES="ceres"       # Skip CERES
CHECKOUT="checkout" # Skip git checkout step
INSTALL=""          # Install as well as build
QUIT_INSTALL_DIR="" # User can specify install directory (mainly for Travis)
NUM_THREADS="2"     # Specify number of threads during make/ninja
NATIVE=""           # Specifieds march=native
while getopts "cij:np:" opt; do
    case $opt in
        c) CERES="";;
        g) CHECKOUT="";;
        i) INSTALL="install";;
        j) NUM_THREADS="$OPTARG";;
        n) NATIVE="-DCMAKE_CXX_FLAGS=-march=native";;
        p) QUIT_INSTALL_DIR="-DCMAKE_INSTALL_PREFIX=$OPTARG";;
    esac
done

# Check for presence of ninja
if [ -x "$(which ninja)" ]; then
    GENERATOR="-GNinja"
    BUILDCMD="ninja -j $NUM_THREADS"
else
    GENERATOR=""
    BUILDCMD="make -j $NUM_THREADS"
fi

# Initialise submodules
if [ -n "$CHECKOUT" ]; then
    git submodule update --init
fi
EXTERNAL="$WD/External"
cd $EXTERNAL

# These are header only, no building required
EIGEN_DIR="$EXTERNAL/eigen"
ARGS_DIR="$EXTERNAL/args"
BUILD_DIR="_build"
CXX_STANDARD="-DCMAKE_CXX_STANDARD=17"

# Ceres
CERES_DIR="$EXTERNAL/ceres-solver"
CERES_BUILD_DIR="$CERES_DIR/$BUILD_DIR"
CERES_OPTS="$GENERATOR -DCMAKE_BUILD_TYPE=Release $NATIVE $CXX_STANDARD\
            -DBUILD_DOCUMENTATION=OFF -DBUILD_EXAMPLES=OFF\
            -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF\
            -DCUSTOM_BLAS=OFF -DCXSPARSE=OFF\
            -DEIGEN_INCLUDE_DIR_HINTS=$EIGEN_DIR -DEIGENSPARSE=OFF\
            -DEXPORT_BUILD_DIR=ON\
            -DGFLAGS=OFF -DLAPACK=OFF -DMINIGLOG=ON\
            -DSUITESPARSE=OFF -DACCELERATESPARSE=OFF"
if [ -n "$CERES" ]; then
    mkdir -p $CERES_BUILD_DIR
    cd $CERES_BUILD_DIR
    cmake $CERES_DIR $CERES_OPTS
    $BUILDCMD
fi

# ITK
ITK_DIR="$EXTERNAL/ITK"
ITK_BUILD_DIR="$ITK_DIR/$BUILD_DIR"
ITK_OPTS="$GENERATOR -DCMAKE_BUILD_TYPE=Release $NATIVE -DCMAKE_CXX_STANDARD=14 -DCMAKE_CXX_FLAGS=-fpermissive
           -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF -DGDCM_USE_COREFOUNDATION_LIBRARY=OFF \
           -DITK_BUILD_DEFAULT_MODULES=OFF\
           -DITKGroup_Core=OFF -DModule_ITKCommon=ON\
           -DModule_ITKIONIFTI=ON -DModule_ITKIONRRD=ON\
           -DModule_ITKIOTransformInsightLegacy=ON\
           -DModule_ITKBinaryMathematicalMorphology=ON\
           -DModule_ITKConnectedComponents=ON\
           -DModule_ITKFFT=ON\
           -DModule_ITKIOImageBase=ON\
           -DModule_ITKIOTransformBase=ON\
           -DModule_ITKImageCompose=ON\
           -DModule_ITKImageFeature=ON\
           -DModule_ITKImageFilterBase=ON\
           -DModule_ITKImageFunction=ON\
           -DModule_ITKImageGrid=ON\
           -DModule_ITKImageIntensity=ON\
           -DModule_ITKImageStatistics=ON\
           -DModule_ITKLabelMap=ON\
           -DModule_ITKLabelVoting=ON\
           -DModule_ITKMathematicalMorphology=ON\
           -DModule_ITKOptimizers=ON\
           -DModule_ITKRegistrationCommon=ON\
           -DModule_ITKSmoothing=ON\
           -DModule_ITKThresholding=ON\
           -DModule_ITKTransform=ON"
mkdir -p $ITK_BUILD_DIR
cd $ITK_BUILD_DIR
cmake $ITK_DIR $ITK_OPTS
$BUILDCMD

cd $WD
# Now build QUIT

QUIT_BLD_DIR="$BUILD_DIR"
BUILDCMD="$BUILDCMD $INSTALL"
QUIT_OPTS="$GENERATOR -DCMAKE_BUILD_TYPE=Release $NATIVE\
    -DEIGEN3_INCLUDE_DIR=$EIGEN_DIR\
    -DCeres_DIR=$CERES_BUILD_DIR\
    -DITK_DIR=$ITK_BUILD_DIR\
    $CXX_STANDARD\
    $QUIT_INSTALL_DIR"
mkdir -p $QUIT_BLD_DIR
cd $QUIT_BLD_DIR
cmake $WD $QUIT_OPTS
$BUILDCMD
