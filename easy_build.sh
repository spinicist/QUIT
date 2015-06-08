#!/bin/bash -ex

#
# easy_build.sh
#
# A hopefully simple script to download all dependencies for QUIT and then
# compile automatically, for those who are unfamiliar with CMake etc.
#
# Requires:
#
# 1. CMake version 2.8 or greater
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

# ITK
ITK_MAJOR="4.7"
ITK_MINOR="2"
ITK_VER="${ITK_MAJOR}.${ITK_MINOR}"
ITK_SRC_DIR="ITK-${ITK_VER}"
ITK_BLD_DIR="ITK-${ITK_VER}-Build"
ITK_URL="http://downloads.sourceforge.net/project/itk/itk/4.7/InsightToolkit-${ITK_VER}.tar.gz?r=&ts=1433761469&use_mirror=kent"
ITK_URL="http://downloads.sourceforge.net/project/itk/itk/4.7/InsightToolkit-${ITK_VER}.tar.gz"
ITK_OPTS="-DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF"
curl --location ${ITK_URL} > ${ITK_DIR}.tar.gz
mkdir -p $ITK_SRC_DIR
tar --extract --file=${ITK_DIR}.tar.gz --strip-components=1 --directory=${ITK_SRC_DIR}
mkdir -p $ITK_BLD_DIR
cd $ITK_BLD_DIR
cmake ../${ITK_SRC_DIR} ${ITK_OPTS}
make -j 2

cd $WD
# Now build QUIT
QUIT_BLD_DIR="Build"
QUIT_OPTS="-DITK_DIR=${WD}/${EXT_DIR}/${ITK_BLD_DIR} -DEIGEN3_INCLUDE_DIR=${WD}/${EXT_DIR}/${EIGEN_DIR}"
mkdir -p $QUIT_BLD_DIR
cd $QUIT_BLD_DIR
cmake $WD ${QUIT_OPTS}
make -j
