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
# First download Eigen & ITK
EXT_DIR="External"
mkdir -p $EXT_DIR
cd $EXT_DIR

# Eigen
EIGEN_VER="3.2.4"
EIGEN_DIR="eigen${EIGEN_VER}"
EIGEN_URL="http://bitbucket.org/eigen/eigen/get/${EIGEN_VER}.tar.gz"
curl --location $EIGEN_URL > ${EIGEN_DIR}.tar.gz
mkdir -p $EIGEN_DIR
tar --extract --file=${EIGEN_DIR}.tar.gz --strip-components=1 --directory=${EIGEN_DIR}
rm ${EIGEN_DIR}.tar.gz

# CppSolvers
SOLVERS_DIR="CppNumericalSolvers"
SOLVERS_URL="http://github.com/spinicist/CppNumericalSolvers/archive/cmaes_rebase.tar.gz"
curl --location $SOLVERS_URL > ${SOLVERS_DIR}.tar.gz
mkdir -p ${SOLVERS_DIR}
tar --extract --file=${SOLVERS_DIR}.tar.gz --strip-components=1 --directory=${SOLVERS_DIR}
rm ${SOLVERS_DIR}.tar.gz

# ITK
ITK_MAJOR="4.10"
ITK_MINOR="0"
ITK_VER="${ITK_MAJOR}.${ITK_MINOR}"
ITK_SRC_DIR="ITK-${ITK_VER}"
ITK_BLD_DIR="ITK-${ITK_VER}-Build"
ITK_URL="https://github.com/InsightSoftwareConsortium/ITK/archive/v${ITK_MAJOR}.${ITK_MINOR}.tar.gz"
ITK_OPTS="-DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF -DCMAKE_CXX_STANDARD=11 -DCMAKE_CXX_FLAGS=-fpermissive"
curl --location ${ITK_URL} > ${ITK_DIR}.tar.gz
mkdir -p $ITK_SRC_DIR
tar --extract --file=${ITK_DIR}.tar.gz --strip-components=1 --directory=${ITK_SRC_DIR}
rm ${ITK_DIR}.tar.gz
mkdir -p $ITK_BLD_DIR
cd $ITK_BLD_DIR
cmake ../${ITK_SRC_DIR} ${ITK_OPTS}
make -j 4

cd $WD
# Now build QUIT
QUIT_BLD_DIR="Build"
QUIT_OPTS="-DCMAKE_BUILD_TYPE=Release -DITK_DIR=${WD}/${EXT_DIR}/${ITK_BLD_DIR} -DEIGEN3_INCLUDE_DIR=${WD}/${EXT_DIR}/${EIGEN_DIR}"
mkdir -p $QUIT_BLD_DIR
cd $QUIT_BLD_DIR
cmake $WD ${QUIT_OPTS}
make -j