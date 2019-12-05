#!/bin/bash -eu

#
# build.sh
#
# A hopefully simple script to download all dependencies for QUIT and then
# compile automatically, for those who are unfamiliar with CMake etc.
#
# Requires:
#
# 1. CMake version 3.2 or greater
# 2. A C++11 compliant compiler e.g. GCC 4.8 or higher
#

USAGE="Usage: $0 [opts]

This script automates the build process for QUIT.

Options:
    -h        Print this help message
    -i        Install after building (default location /usr/local/bin)
    -j N      Specify number of parallel build processes
    -p PREFIX Change the install location prefix (/bin will be added)
    -n        Add -march=native to build flags (improves performance)
    -s gci    Skip (g)it submodule update/(c)eres build/(i)TK build steps
"

WD=$PWD

# Command line options
CHECKOUT="true"     # Run git checkout step
CERES="true"        # Build Ceres
ITK="true"          # Build ITK
INSTALL=""          # Install as well as build
QUIT_INSTALL_DIR="" # User can specify install directory (mainly for Travis)
NUM_THREADS=""      # Specify number of threads during make/ninja
NATIVE=""           # Specifieds march=native
while getopts "hij:np:s:" opt; do
    case $opt in
        h) echo "$USAGE"; exit 0;;
        i) INSTALL="install";;
        j) NUM_THREADS="$OPTARG";;
        n) NATIVE="-DCMAKE_CXX_FLAGS=-march=native";;
        p) QUIT_INSTALL_DIR="-DCMAKE_INSTALL_PREFIX=$OPTARG";;
        s)
            if [[ $OPTARG == *"g"* ]]; then
                CHECKOUT=""
            fi
            if [[ $OPTARG == *"c"* ]]; then
                CERES=""
            fi
            if [[ $OPTARG == *"i"* ]]; then
                ITK=""
            fi
    esac
done
shift $(( $OPTIND - 1))

if [ $# -ne 0 ]; then
    echo "$USAGE"
    exit 1;
fi

# Check for presence of ninja
if [ -x "$(which ninja)" ]; then
    GENERATOR="-GNinja"
    BUILDCMD="ninja -j $NUM_THREADS"
    echo "Found ninja"
else
    GENERATOR=""
    BUILDCMD="make -j $NUM_THREADS"
fi

# Initialise submodules
if [ -n "$CHECKOUT" ]; then
    echo "Running git submodule update"
    git submodule update --init
fi
EXTERNAL="$WD/External"
cd $EXTERNAL

# These are header only, no building required
EIGEN_DIR="$EXTERNAL/eigen"
ARGS_DIR="$EXTERNAL/args"
BUILD_DIR="build"
CXX_STANDARD="-DCMAKE_CXX_STANDARD=17"

# Ceres
CERES_DIR="$EXTERNAL/ceres-solver"
# Needs an underscore because of Bazel BUILD file
CERES_BUILD_DIR="${CERES_DIR}/_${BUILD_DIR}"
if [ -n "$CERES" ]; then
    echo "Building ceres..."
    if [[ -d $CERES_BUILD_DIR ]]; then
        rm -r $CERES_BUILD_DIR
    fi
    mkdir -p "$CERES_BUILD_DIR"
    cd "$CERES_BUILD_DIR"
    CERES_OPTS="$GENERATOR -DCMAKE_BUILD_TYPE=Release $NATIVE $CXX_STANDARD\
    -DBUILD_DOCUMENTATION=OFF -DBUILD_EXAMPLES=OFF\
    -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF\
    -DCUSTOM_BLAS=OFF -DCXSPARSE=OFF\
    -DEIGEN_INCLUDE_DIR_HINTS=$EIGEN_DIR -DEIGENSPARSE=OFF\
    -DEXPORT_BUILD_DIR=ON\
    -DGFLAGS=OFF -DLAPACK=OFF -DMINIGLOG=ON -DMINIGLOG_MAX_LOG_LEVEL=-100\
    -DSUITESPARSE=OFF -DACCELERATESPARSE=OFF"
    cmake $CERES_DIR $CERES_OPTS
    $BUILDCMD
fi

# ITK
ITK_DIR="$EXTERNAL/ITK"
ITK_BUILD_DIR="$ITK_DIR/$BUILD_DIR"
if [ -n "$ITK" ]; then
    echo "Building ITK..."
    if [[ -d $ITK_BUILD_DIR ]]; then
        rm -r $ITK_BUILD_DIR
    fi
    mkdir -p "$ITK_BUILD_DIR"
    cd "$ITK_BUILD_DIR"
    ITK_OPTS="$GENERATOR -DCMAKE_BUILD_TYPE=Release $NATIVE -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_FLAGS=-fpermissive
    -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF \
    -DITK_BUILD_DEFAULT_MODULES=OFF\
    -DITKGroup_Core=OFF\
    -DModule_ITKBinaryMathematicalMorphology=ON\
    -DModule_ITKCommon=ON\
    -DModule_ITKConnectedComponents=ON\
    -DModule_ITKFFT=ON\
    -DModule_ITKIOImageBase=ON\
    -DModule_ITKIOTransformBase=ON\
    -DModule_ITKImageCompose=ON\
    -DModule_ITKImageFeature=ON\
    -DModule_ITKImageFilterBase=ON\
    -DModule_ITKImageGrid=ON\
    -DModule_ITKImageIntensity=ON\
    -DModule_ITKImageStatistics=ON\
    -DModule_ITKLabelMap=ON\
    -DModule_ITKLabelVoting=ON\
    -DModule_ITKMathematicalMorphology=ON\
    -DModule_ITKThresholding=ON\
    -DModule_ITKTransform=ON\
    -DModule_ITKIOTransformInsightLegacy=ON\
    -DModule_ITKIONIFTI=ON"
    cmake $ITK_DIR $ITK_OPTS
    $BUILDCMD
fi

# Now build QUIT
QUIT_BUILD_DIR="$BUILD_DIR"
BUILDCMD="$BUILDCMD $INSTALL"
QUIT_OPTS="$GENERATOR -DCMAKE_BUILD_TYPE=Release $NATIVE\
    -DEIGEN3_INCLUDE_DIR=$EIGEN_DIR\
    -DCERES_BUILD_DIR="$CERES_BUILD_DIR"
    -DITK_DIR=$ITK_BUILD_DIR\
    $CXX_STANDARD\
    $QUIT_INSTALL_DIR"
cd $WD
if [[ -d "$QUIT_BUILD_DIR" ]]; then
    rm -r $QUIT_BUILD_DIR
fi
mkdir -p $QUIT_BUILD_DIR
cd $QUIT_BUILD_DIR
echo "Building QUIT..."
cmake $WD $QUIT_OPTS
$BUILDCMD
