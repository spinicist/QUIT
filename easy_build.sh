#!/bin/bash -ex

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

# First initialize submodules
EXT_DIR="$WD/External"
cd $EXT_DIR
git submodule init && git submodule update

# Eigen is header only

# ITK
ITK_DIR="$EXTERNAL/ITK"
ITK_BUILD_DIR="$ITK_DIR/build"
ITK_OPTS="-DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF\
 -DModule_ITKReview=On -DCMAKE_CXX_STANDARD=11 -DCMAKE_CXX_FLAGS=-fpermissive"
mkdir -p $ITK_BUILD_DIR
cd $ITK_BUILD_DIR
cmake $ITK_DIR $ITK_OPTS
make -j 4

# Ceres
CERES_DIR="$EXTERNAL/ceres-solver"
CERES_BUILD_DIR="$CERES_DIR/build"
CERES_OPTS="-DCMAKE_BUILD_TYPE=Release\
            -DBUILD_DOCUMENTATION=OFF -DBUILD_EXAMPLES=OFF\
            -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF\
            -DCUSTOM_BLAS=OFF -DCXSPARSE=OFF -DCXX11=ON\
            -DEIGENSPARSE=OFF -DEXPORT_BUILD_DIR=ON\
            -DGFLAGS=OFF -DLAPACK=OFF -DMINIGLOG=ON\
            -DOPENMP=OFF -DSUITESPARSE=OFF"

mkdir -p $CERES_BUILD_DIR
cd $CERES_BUILD_DIR
cmake $CERES_DIR $CERES_OPTS

cd $WD
# Now build QUIT
QUIT_BLD_DIR="$WD/build"
QUIT_OPTS="-DCMAKE_BUILD_TYPE=Release -DITK_DIR=${WD}/${EXT_DIR}/${ITK_BLD_DIR} -DEIGEN3_INCLUDE_DIR=${WD}/${EXT_DIR}/${EIGEN_DIR}"
mkdir -p $QUIT_BLD_DIR
cd $QUIT_BLD_DIR
cmake $WD ${QUIT_OPTS}
make -j 4