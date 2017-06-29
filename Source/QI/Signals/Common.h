/*
 *  Common.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SIGNALS_COMMON_H
#define SIGNALS_COMMON_H

#include <iostream>
#include <exception>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace QI {

//******************************************************************************
// Magnetisation Evolution Matrices, helper functions etc.
//******************************************************************************
typedef const double cdbl; // To save tedious typing
typedef const Eigen::ArrayXd carrd;

typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 9, 9> Matrix9d;
typedef Eigen::Matrix<double, 9, 1> Vector9d;
typedef Eigen::Matrix<double, 3, Eigen::Dynamic> MagVector;

Eigen::VectorXd SigMag(const MagVector &M_in);
Eigen::VectorXcd SigComplex(const MagVector &M_in);
MagVector SumMC(const Eigen::MatrixXd &M_in);

const Eigen::Matrix3d Relax(cdbl &T1, cdbl &T2);
const Eigen::Matrix3d InfinitesimalRF(cdbl &dalpha);
const Eigen::Matrix3d OffResonance(cdbl &inHertz);
const Eigen::Matrix3d Spoiling();
const Matrix6d Exchange(cdbl &k_ab, cdbl &k_ba);
const void CalcExchange(cdbl tau_a, cdbl f_a, double &f_b, double &k_ab, double &k_ba);

} // End namespace QI

#endif // SIGNALS_COMMON_H