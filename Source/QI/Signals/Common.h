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
#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace std;
using namespace Eigen;

//******************************************************************************
// Convenience
//******************************************************************************
double clamp(double value, double low, double high);

//******************************************************************************
// Magnetisation Evolution Matrices, helper functions etc.
//******************************************************************************
typedef const double cdbl; // To save tedious typing
typedef const ArrayXd carrd;

typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 9, 9> Matrix9d;
typedef Matrix<double, 9, 1> Vector9d;
typedef Matrix<double, 3, Dynamic> MagVector;

VectorXd SigMag(const MagVector &M_in);
VectorXcd SigComplex(const MagVector &M_in);
MagVector SumMC(const MatrixXd &M_in);

const Matrix3d Relax(cdbl &T1, cdbl &T2);
const Matrix3d InfinitesimalRF(cdbl &dalpha);
const Matrix3d OffResonance(cdbl &inHertz);
const Matrix3d Spoiling();
const Matrix6d Exchange(cdbl &k_ab, cdbl &k_ba);
const void CalcExchange(cdbl tau_a, cdbl f_a, double &f_b, double &k_ab, double &k_ba);

#endif // SIGNALS_COMMON_H