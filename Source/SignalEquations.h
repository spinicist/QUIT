/*
 *  SignalEquations.h
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright (c) 2011-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef DESPOT_SIGEQU
#define DESPOT_SIGEQU

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

inline const Matrix3d Relax(cdbl &T1, cdbl &T2);
inline const Matrix3d InfinitesimalRF(cdbl &dalpha);
inline const Matrix3d OffResonance(cdbl &inHertz);
inline const Matrix3d Spoiling();
inline const Matrix6d Exchange(cdbl &k_ab, cdbl &k_ba);
const void CalcExchange(cdbl tau_a, cdbl f_a, double &f_b, double &k_ab, double &k_ba);

//******************************************************************************
// Actual Signal Equations
//******************************************************************************
VectorXcd One_MultiEcho(carrd &TE, cdbl PD, cdbl T1);

VectorXcd One_SPGR(carrd &flip, cdbl TR, cdbl PD, cdbl T1, cdbl B1);
VectorXcd One_SPGR_Echo(carrd &flip, cdbl TR, cdbl TE, cdbl PD, cdbl T1, cdbl T2, cdbl B1);
VectorXcd One_SSFP(carrd &flip, cdbl TR, cdbl ph, cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1);
VectorXcd One_SSFP_Echo(carrd &flip, cdbl TR, cdbl ph, cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1);
VectorXcd One_SSFP_Finite(carrd &flip, const bool spoil, cdbl TR, cdbl Trf, cdbl TE, cdbl ph,
                          cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1);
VectorXcd One_SSFP_Ellipse(carrd &flip, cdbl TR, cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1);
VectorXcd MP_RAGE(cdbl flip, cdbl TR, const int N, carrd &TI, cdbl TD, cdbl PD, cdbl T1, cdbl B1);
VectorXcd One_AFI(cdbl flip, cdbl TR1, cdbl TR2, cdbl PD, cdbl T1, cdbl B1);

VectorXcd Two_SPGR(carrd &flip, cdbl TR, cdbl PD, cdbl T1_a, cdbl T1_b, cdbl tau_a, cdbl f_a, cdbl B1);
VectorXcd Two_SPGR_Echo(carrd &flip, cdbl TR, cdbl TE, cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b, cdbl tau_a, cdbl f_a, cdbl B1);
VectorXcd Two_SSFP(carrd &flip, cdbl TR, cdbl ph, cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b, cdbl tau_a, cdbl f_a, cdbl f0_a, cdbl f0_b, cdbl B1);
VectorXcd Two_SSFP_Echo(carrd &flip, cdbl TR, cdbl ph, cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b, cdbl tau_a, cdbl f_a, cdbl f0_a, cdbl f0_b, cdbl B1);
VectorXcd Two_SSFP_Finite(carrd &flip, const bool spoil, cdbl TR, cdbl Trf, cdbl TE, cdbl ph,
                          cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b,
                          cdbl tau_a, cdbl f_a, cdbl f0_a, cdbl f0_b, cdbl B1);

VectorXcd Three_SPGR(carrd &flip, cdbl TR, cdbl PD, cdbl T1_a, cdbl T1_b, cdbl T1_c, cdbl tau_a, cdbl f_a, cdbl f_c, cdbl B1);
VectorXcd Three_SSFP(carrd &flip, cdbl TR, cdbl ph, cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b, cdbl T1_c, cdbl T2_c, cdbl tau_a, cdbl f_a, cdbl f_c, cdbl f0_a, cdbl f0_b, cdbl f0_c, cdbl B1);
VectorXcd Three_SSFP_Echo(carrd &flip, cdbl TR, cdbl ph, cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b, cdbl T1_c, cdbl T2_c, cdbl tau_a, cdbl f_a, cdbl f_c, cdbl f0_a, cdbl f0_b, cdbl f0_c, cdbl B1);
VectorXcd Three_SSFP_Finite(carrd &flip, const bool spoil, cdbl TR, cdbl Trf, cdbl TE, cdbl ph,
                            cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b, cdbl T1_c, cdbl T2_c,
                            cdbl tau_a, cdbl f_a, cdbl f_c,
                            cdbl f0_a, cdbl f0_b, cdbl f0_c, cdbl B1);

#endif
