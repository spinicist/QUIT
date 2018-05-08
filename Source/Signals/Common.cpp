/*
 *  Common.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Common.h"

using namespace std;
using namespace Eigen;

namespace QI {

//******************************************************************************
#pragma mark Magnetisation Evolution Matrices, helper functions etc.
//******************************************************************************
typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 9, 9> Matrix9d;
typedef Matrix<double, 9, 1> Vector9d;
typedef Matrix<double, 3, Dynamic> MagVector;

// Sum a multi-component magnetisation vector
MagVector SumMC(const MatrixXd &M_in) {
    MagVector M_out = M_in.topRows(3);
    for (MatrixXd::Index i = 3; i < M_in.rows(); i += 3) {
        M_out += M_in.block(i, 0, 3, M_in.cols());
    }
    return M_out;
}

// Turn a full XYZ magnetisation vector into a complex transverse magnetisation
VectorXcd SigComplex(const MagVector &M_in) {
    VectorXcd cs(M_in.cols());
    cs.real() = M_in.topRows(1).transpose();
    cs.imag() = M_in.block(1, 0, 1, M_in.cols()).transpose();
    return cs;
}

VectorXd SigMag(const MagVector &M_in) {
    VectorXd s = M_in.topRows(2).colwise().norm();
    return s;
}

const Matrix3d Relax(const double &T1, const double &T2) {
    Matrix3d R;
    R << 1./T2,     0,     0,
             0, 1./T2,     0,
             0,     0, 1./T1;
    return R;
}

const Matrix3d InfinitesimalRF(const double &dalpha) {
    Matrix3d A;
    A << 0,      0, -dalpha,
         0,      0,       0,
         dalpha, 0,       0;
    return A;
}

const Matrix3d OffResonance(const double &inHertz) {
    // Minus signs are this way round to make phase cycling go the right way
    Matrix3d O;
    double dw = inHertz * 2. * M_PI;
    O <<  0, dw, 0,
        -dw,  0, 0,
          0,  0, 0;
    return O;
}

const Matrix3d Spoiling() {
    // Destroy the x- and y- magnetization
    Matrix3d S = Matrix3d::Zero();
    S(2, 2) = 1.;
    return S;
}

const Matrix6d Exchange(const double &k_ab, const double &k_ba) {
    Matrix6d K = Matrix6d::Zero();
    K.block(0,0,3,3).diagonal().setConstant(k_ab);
    K.block(3,3,3,3).diagonal().setConstant(k_ba);
    K.block(3,0,3,3).diagonal().setConstant(-k_ab);
    K.block(0,3,3,3).diagonal().setConstant(-k_ba);
    return K;
}

// Calculate the exchange rates from the residence time and fractions
void CalcExchange(const double tau_a, const double f_a, double &f_b, double &k_ab, double &k_ba) {
    const double feps = numeric_limits<float>::epsilon(); // Because we read from float files
    f_b = 1.0 - f_a;
    k_ab = 1./tau_a; k_ba = k_ab*f_a/f_b;
    if ((fabs(f_a - 1.) <= feps) || (fabs(f_b - 1.) <= feps)) {
        // Only have 1 component, so no exchange
        k_ab = 0.;
        k_ba = 0.;
    }
}

} // End namespace QI
