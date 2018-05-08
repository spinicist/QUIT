/*
 *  SPGR.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "SPGR.h"
// There's a bug in matrix log that leads to duplicate symbol definitions if
// this is included in the header file
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;
using namespace Eigen;

namespace QI {

VectorXcd One_SPGR(carrd &flip, cdbl TR, cdbl PD, cdbl T1, cdbl B1) {
    VectorXcd M = VectorXcd::Zero(flip.size());
    ArrayXd sa = (B1 * flip).sin();
    ArrayXd ca = (B1 * flip).cos();
    double expT1 = exp(-TR / T1);
    M.real() = PD * ((1. - expT1) * sa) / (1. - expT1*ca);
    return M;
}

VectorXd One_SPGR_Magnitude(carrd &flip, cdbl TR, cdbl PD, cdbl T1, cdbl B1) {
    const ArrayXd sa = sin(B1 * flip);
    const ArrayXd ca = cos(B1 * flip);
    const double E1 = exp(-TR / T1);
    return PD * ((1. - E1) * sa) / (1. - E1*ca);
}

ArrayXXd One_SPGR_Magnitude_Derivs(carrd &flip, cdbl TR, cdbl PD, cdbl T1, cdbl B1) {
    const ArrayXd sa = sin(B1 * flip);
    const ArrayXd ca = cos(B1 * flip);
    const double E1 = exp(-TR / T1);
    const ArrayXd d = (1.-E1*ca);
    ArrayXXd derivs(flip.rows(), 2);
    derivs.col(0) = (1.-E1)*sa/d;
    derivs.col(1) = E1*PD*TR*(ca-1.)*sa/((d*T1).square());
    return derivs;
}


VectorXcd One_SPGR_Echo(carrd &flip, cdbl TR, cdbl TE, cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1) {
    ArrayXd sa = (B1 * flip).sin();
    ArrayXd ca = (B1 * flip).cos();
    double E1 = exp(-TR / T1);
    double E2 = exp(-TE / T2); // T2' is incorporated into PD in this formalism
    complex<double> phase = polar(PD * E2, 2. * M_PI * TE * f0);
    VectorXcd M = phase * (((1. - E1) * sa) / (1. - E1*ca));
    return M;
}

VectorXcd Two_SPGR(carrd &flip, cdbl TR,
                   cdbl PD, cdbl T1_a, cdbl T1_b, cdbl tau_a, cdbl f_a, cdbl B1) {
    Matrix2d A, eATR;
    Vector2d M0, Mobs;
    MagVector signal(3, flip.size()); signal.setZero();
    double k_ab, k_ba, f_b;
    CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    M0 << f_a, f_b;
    A << ((1./T1_a) + k_ab),            -k_ba,
                     -k_ab, ((1./T1_b) + k_ba);
    eATR = (-TR*A).exp();
    const Vector2d RHS = (Matrix2d::Identity() - eATR) * M0;
    for (int i = 0; i < flip.size(); i++) {
        const double a = flip[i] * B1;
        Mobs = (Matrix2d::Identity() - eATR*cos(a)).partialPivLu().solve(RHS * sin(a));
        signal(1, i) = PD * Mobs.sum();
    }
    return SigComplex(signal);
}

VectorXcd Two_SPGR_Echo(carrd &flip, cdbl TR, cdbl TE, cdbl PD,
                        cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b, 
                        cdbl tau_a, cdbl f_a, 
                        cdbl f0_a, cdbl f0_b, cdbl B1) {
    Matrix2d A, eATR;
    Vector2d M0, Mz;
    Vector2cd Mxy;
    double k_ab, k_ba, f_b;
    CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    M0 << f_a, f_b;
    A << -(1./T1_a) - k_ab,              k_ba,
                      k_ab, -(1./T1_b) - k_ba;
    eATR = (TR*A).exp();
    Matrix2cd echo = Matrix2cd::Zero();
    echo(0,0) = polar(exp(-TE/T2_a),2.*M_PI*f0_a); // T2' absorbed into PD as it effects both components equally
    echo(1,1) = polar(exp(-TE/T2_b),2.*M_PI*f0_b);
    const Vector2d RHS = (Matrix2d::Identity() - eATR) * M0;
    VectorXcd signal(flip.size());
    for (int i = 0; i < flip.size(); i++) {
        const double a = flip[i] * B1;
        Mz.noalias() = (Matrix2d::Identity() - eATR*cos(a)).partialPivLu().solve(RHS);
        Mxy.noalias() = echo * Mz * sin(a);
        signal(i) = PD * (Mxy(0) + Mxy(1));
    }
    return signal;
}

VectorXcd Three_SPGR(carrd &flip, cdbl TR, cdbl PD,
                     cdbl T1_a, cdbl T1_b, cdbl T1_c, cdbl tau_a, cdbl f_a, cdbl f_c, cdbl B1) {
    double f_ab = 1. - f_c;
    VectorXcd m_ab = Two_SPGR(flip, TR, PD * f_ab, T1_a, T1_b, tau_a, f_a / f_ab, B1);
    VectorXcd m_c  = One_SPGR(flip, TR, PD * f_c, T1_c, B1);
    VectorXcd r = m_ab + m_c;
    return r;
}

VectorXcd Three_SPGR_Echo(carrd &flip, cdbl TR, cdbl TE, cdbl PD,
                          cdbl T1_a, cdbl T2_a, 
                          cdbl T1_b, cdbl T2_b, 
                          cdbl T1_c, cdbl T2_c, 
                          cdbl tau_a, cdbl f_a, cdbl f_c, 
                          cdbl f0_a, cdbl f0_b, cdbl f0_c,
                          cdbl B1) {
    double f_ab = 1. - f_c;
    VectorXcd m_ab = Two_SPGR_Echo(flip, TR, TE, PD * f_ab, T1_a, T2_a, T1_b, T2_b, tau_a, f_a / f_ab, f0_a, f0_b, B1);
    VectorXcd m_c  = One_SPGR_Echo(flip, TR, TE, PD * f_c, T1_c, T2_c, f0_c, B1);
    VectorXcd r = m_ab + m_c;
    return r;
}

VectorXcd MT_SPGR(carrd &omega_cwpe, carrd &satf0, cdbl /* Unused */, cdbl /* Unused */,
                  const TLineshape &g, cdbl PD, cdbl T1f, cdbl T2f, cdbl T1r, cdbl T2r,
                  cdbl kf, cdbl F, cdbl /* Unused */, cdbl /* Unused */) {
    // This is a mix of Gloor et al and Ramani et al variable names
    ArrayXd W = omega_cwpe.square()*g(satf0, T2r);
    cdbl R1f = 1. / T1f;
    cdbl R1r = 1. / T1r;
    
    // F is M0r/M0b
    cdbl kr = kf/F;
    carrd S = PD * F * ( R1r*kr/R1f + W + R1r + kr ) /
                    ( kf*(R1r + W) + (1.0 + (omega_cwpe/(2.*M_PI*satf0)).square()*(T1f/T2f))*(W+R1r+kr));
    VectorXcd sig(S.size());
    sig.real() = S;
    return sig;
}

} // End namespace QI
