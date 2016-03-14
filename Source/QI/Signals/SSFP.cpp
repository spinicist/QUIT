/*
 *  SSFP.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Signals/SSFP.h"
// There's a bug in matrix log that leads to duplicate symbol definitions if
// this is included in the header file
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;
using namespace Eigen;

namespace QI {

VectorXcd One_SSFP_GS(carrd &flip, cdbl TR, cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1) {
    // This is at the echo time
    const double E1 = exp(-TR / T1);
    const double E2 = exp(-TR / T2);

    const double psi = 2. * M_PI * f0 * TR;
    const ArrayXd alpha = flip * B1;
    const ArrayXd d = (1. - E1*E2*E2-(E1-E2*E2)*cos(alpha));
    const ArrayXcd G = -PD*sqrt(E2)*polar(1., psi/2)* (1 - E1)*sin(alpha)/d;
    return G;
}

VectorXcd One_SSFP(carrd &flip, carrd &phi, cdbl TR,
                   cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1) {
    eigen_assert(flip.size() == phi.size())
    const double E1 = exp(-TR / T1);
    const double E2 = exp(-TR / T2);

    const double  psi = 2. * M_PI * f0 * TR;
    const ArrayXd alpha = flip * B1;
    const ArrayXd theta = phi + psi;
    const ArrayXd d = (1. - E1*E2*E2-(E1-E2*E2)*cos(alpha));
    const ArrayXd G = -PD*(1. - E1)*sin(alpha)/d;
    const ArrayXd b = E2*(1. - E1)*(1.+cos(alpha))/d;
    ArrayXcd et(theta.size());
    et.real() = cos(-theta);
    et.imag() = sin(-theta);
    const ArrayXcd M = G*(1. - E2*et) / (1 - b*cos(theta));
    return M;
}

VectorXcd One_SSFP_Echo(carrd &flip, carrd &phi, cdbl TR,
                        cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1) {
    eigen_assert(flip.size() = phi.size());
    const double E1 = exp(-TR / T1);
    const double E2 = exp(-TR / T2);

    const double   psi = 2. * M_PI * f0 * TR;
    const ArrayXd  alpha = flip * B1;
    const ArrayXd  theta = phi + psi;
    const ArrayXd  d = (1. - E1*E2*E2-(E1-E2*E2)*cos(alpha));
    const ArrayXcd G = -PD*sqrt(E2)*polar(1., psi/2.)*(1 - E1)*sin(alpha)/d;
    const ArrayXd  b = E2*(1. - E1)*(1.+cos(alpha))/d;
    ArrayXcd et(theta.size());
    et.real() = cos(-theta);
    et.imag() = sin(-theta);
    const ArrayXcd M = G*(1. - E2*et) / (1 - b*cos(theta));
    return M;
}

VectorXcd One_SSFP_Finite(carrd &flip, const bool spoil, cdbl TR, cdbl Trf, cdbl inTE, cdbl phase,
                          cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1) {
    const Matrix3d I = Matrix3d::Identity();
    const Matrix3d O = OffResonance(f0);
    Matrix3d P, R = Relax(T1, T2);
    double TE;
    if (spoil) {
        P = Spoiling();
        TE = inTE - Trf;
        assert(TE > 0.);
    } else {
        P = AngleAxisd(phase, Vector3d::UnitZ());
        TE = (TR - Trf) / 2; // Time AFTER the RF Pulse ends that echo is formed
    }
    
    const Matrix3d RpO = R + O;
    const Matrix3d E_e = (-TE * RpO).exp();
    const Matrix3d E = (-(TR - Trf) * RpO).exp();
    Vector3d m_inf; m_inf << 0, 0, PD;
        
    Matrix3d E_r;
    MagVector result(3, flip.size());
    for (int i = 0; i < flip.size(); i++) {
        const Matrix3d A = InfinitesimalRF(B1 * flip(i) / Trf);
        E_r.noalias() = (-Trf * (RpO + A)).exp();
        Vector3d m_rinf = (RpO + A).partialPivLu().solve(R * m_inf);
        Vector3d m_r = (I - E_r*P*E).partialPivLu().solve(E_r*P*(I-E)*m_inf + (I-E_r)*m_rinf);
        Vector3d m_e = E_e*(m_r - m_inf) + m_inf;
        result.col(i) = m_e;
    }
    return SigComplex(result);
}

MatrixXd Two_SSFP_Matrix(carrd &flip, carrd &phi, const double TR,
                         cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b,
                         cdbl tau_a, cdbl f_a, cdbl f0_a, cdbl f0_b, cdbl B1) {
    const double E1_a = exp(-TR/T1_a);
    const double E1_b = exp(-TR/T1_b);
    const double E2_a = exp(-TR/T2_a);
    const double E2_b = exp(-TR/T2_b);
    double f_b, k_ab, k_ba;
    CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    const double E_ab = exp(-TR*k_ab/f_b);
    const double K1 = E_ab*f_b+f_a;
    const double K2 = E_ab*f_a+f_b;
    const double K3 = f_a*(1-E_ab);
    const double K4 = f_b*(1-E_ab);
    carrd alpha = B1 * flip;
    carrd theta_a = phi + 2.*M_PI*f0_a*TR;
    carrd theta_b = phi + 2.*M_PI*f0_b*TR;

    MatrixXd M(4, flip.size());
    Matrix6d LHS;
    Vector6d RHS;
    RHS << 0, 0, 0, 0, -E1_b*K3*f_b + f_a*(-E1_a*K1 + 1), -E1_a*K4*f_a + f_b*(-E1_b*K2 + 1);

    for (int i = 0; i < flip.size(); i++) {
        const double ca = cos(alpha[i]);
        const double sa = sin(alpha[i]);
        const double cta = cos(theta_a[i]);
        const double ctb = cos(theta_b[i]);
        const double sta = sin(theta_a[i]);
        const double stb = sin(theta_b[i]);

        LHS << -E2_a*K1*cta + ca, -E2_b*K3*cta, E2_a*K1*sta, E2_b*K3*sta, sa, 0,
               -E2_a*K4*ctb, -E2_b*K2*ctb + ca, E2_a*K4*stb, E2_b*K2*stb, 0, sa,
               -E2_a*K1*sta, -E2_b*K3*sta, -E2_a*K1*cta + 1, -E2_b*K3*cta, 0, 0,
               -E2_a*K4*stb, -E2_b*K2*stb, -E2_a*K4*ctb, -E2_b*K2*ctb + 1, 0, 0,
               -sa, 0, 0, 0, -E1_a*K1 + ca, -E1_b*K3,
                0, -sa, 0, 0, -E1_a*K4, -E1_b*K2 + ca;
        M.col(i).noalias() = (LHS.partialPivLu().solve(RHS)).head(4);
    }
    return M;
}

VectorXcd Two_SSFP(carrd &flip, carrd &phi, const double TR,
                   cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b,
                   cdbl tau_a, cdbl f_a, cdbl f0_a, cdbl f0_b, cdbl B1) {
    eigen_assert(flip.size() == phi.size())

    MatrixXd M = Two_SSFP_Matrix(flip, phi, TR, T1_a, T2_a, T1_b, T2_b, tau_a, f_a, f0_a, f0_b, B1);
    VectorXcd mc(flip.size());
    mc.real() = PD * (M.row(0) + M.row(1));
    mc.imag() = PD * (M.row(2) + M.row(3));
    return mc;
}

VectorXcd Two_SSFP_Echo(carrd &flip, carrd &phi, const double TR,
                        cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b,
                        cdbl tau_a, cdbl f_a, cdbl f0_a, cdbl f0_b, cdbl B1) {
    eigen_assert(flip.size() == phi.size())

    MatrixXd M = Two_SSFP_Matrix(flip, phi, TR, T1_a, T2_a, T1_b, T2_b, tau_a, f_a, f0_a, f0_b, B1);

    const double sE2_a = exp(-TR/(2.*T2_a));
    const double sE2_b = exp(-TR/(2.*T2_b));
    double f_b, k_ab, k_ba;
    CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    const double sqrtE_ab = exp(-TR*k_ab/(2*f_b));
    const double K1 = sqrtE_ab*f_b+f_a;
    const double K2 = sqrtE_ab*f_a+f_b;
    const double K3 = f_a*(1-sqrtE_ab);
    const double K4 = f_b*(1-sqrtE_ab);
    const double cpa = cos(M_PI*f0_a*TR);
    const double spa = sin(M_PI*f0_a*TR);
    const double cpb = cos(M_PI*f0_b*TR);
    const double spb = sin(M_PI*f0_b*TR);

    Matrix<double, 4, 4> echo;
    echo << sE2_a*K1*cpa, sE2_b*K3*cpa, -sE2_a*K1*spa, -sE2_b*K3*spa,
            sE2_a*K4*cpb, sE2_b*K2*cpb, -sE2_a*K4*spb, -sE2_b*K2*spb,
            sE2_a*K1*spa, sE2_b*K3*spa,  sE2_a*K1*cpa,  sE2_b*K3*cpa,
            sE2_a*K4*spb, sE2_b*K2*spb,  sE2_a*K4*cpb,  sE2_b*K2*cpb;

    const MatrixXd Me = echo*M;
    VectorXcd mce(flip.size());
    mce.real() = PD * (Me.row(0) + Me.row(1));
    mce.imag() = PD * (Me.row(2) + Me.row(3));
    return mce;
}

VectorXcd Two_SSFP_Finite(carrd &flip, const bool spoil,
                          cdbl TR, cdbl Trf, cdbl inTE, cdbl phase,
                          cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b,
                          cdbl tau_a, cdbl f_a, cdbl f0_a, cdbl f0_b, cdbl B1) {
    const Matrix6d I = Matrix6d::Identity();
    Matrix6d R = Matrix6d::Zero(), C = Matrix6d::Zero();
    Matrix6d O = Matrix6d::Zero();
    O.block(0,0,3,3) = OffResonance(f0_a);
    O.block(3,3,3,3) = OffResonance(f0_b);
    Matrix3d C3;
    double TE;
    if (spoil) {
        TE = inTE - Trf;
        assert(TE > 0.);
        C3 = Spoiling();
        R.block(0,0,3,3) = Relax(T1_a, T2_a);
        R.block(3,3,3,3) = Relax(T1_b, T2_b);
    } else {
        TE = (TR - Trf) / 2; // Time AFTER the RF Pulse ends that echo is formed
        C3 = AngleAxisd(phase, Vector3d::UnitZ());
        R.block(0,0,3,3) = Relax(T1_a, T2_a);
        R.block(3,3,3,3) = Relax(T1_b, T2_b);
    }
    C.block(0,0,3,3) = C.block(3,3,3,3) = C3;
    Matrix6d RpO = R + O;
    double k_ab, k_ba, f_b;
    CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    Matrix6d K = Exchange(k_ab, k_ba);
    Matrix6d RpOpK = RpO + K;
    Matrix6d l1;
    const Matrix6d le = (-(RpOpK)*TE).exp();
    const Matrix6d l2 = (-(RpOpK)*(TR-Trf)).exp();
    
    Vector6d m0, mp, me;
    m0 << 0, 0, f_a * PD, 0, 0, PD * f_b;
    const Vector6d Rm0 = R * m0;
    const Vector6d m2 = (RpO).partialPivLu().solve(Rm0);
    const Vector6d Cm2 = C * m2;
    
    MagVector theory(3, flip.size());
    Matrix6d A = Matrix6d::Zero();
    
    for (int i = 0; i < flip.size(); i++) {
        A.block(0,0,3,3) = A.block(3,3,3,3) = InfinitesimalRF(B1 * flip(i) / Trf);
        l1.noalias() = (-(RpOpK+A)*Trf).exp();
        Vector6d m1 = (RpO + A).partialPivLu().solve(Rm0);
        mp.noalias() = Cm2 + (I - l1*C*l2).partialPivLu().solve((I - l1)*(m1 - Cm2));
        me.noalias() = le*(mp - m2) + m2;				
        theory.col(i) = SumMC(me);
    }
    return SigComplex(theory);
}


VectorXcd Three_SSFP(carrd &flip, carrd &phi, cdbl TR, cdbl PD,
                     cdbl T1_a, cdbl T2_a,
                     cdbl T1_b, cdbl T2_b,
                     cdbl T1_c, cdbl T2_c,
                     cdbl tau_a, cdbl f_a, cdbl f_c,
                     cdbl f0_a, cdbl f0_b, cdbl f0_c, cdbl B1) {
    double f_ab = 1. - f_c;
    VectorXcd m_ab = Two_SSFP(flip, phi, TR, PD * f_ab, T1_a, T2_a, T1_b, T2_b, tau_a, f_a / f_ab, f0_a, f0_b, B1);
    VectorXcd m_c  = One_SSFP(flip, phi, TR, PD * f_c, T1_c, T2_c, f0_c, B1);
    VectorXcd r = m_ab + m_c;
    return r;
}

VectorXcd Three_SSFP_Echo(carrd &flip, carrd &phi, cdbl TR, 
                          cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b, cdbl T1_c, cdbl T2_c,
                          cdbl tau_a, cdbl f_a, cdbl f_c, cdbl f0_a, cdbl f0_b, cdbl f0_c, cdbl B1) {
    double f_ab = 1. - f_c;
    VectorXcd m_ab = Two_SSFP_Echo(flip, phi, TR, PD * f_ab, T1_a, T2_a, T1_b, T2_b, tau_a, f_a / f_ab, f0_a, f0_b, B1);
    VectorXcd m_c  = One_SSFP_Echo(flip, phi, TR, PD * f_c, T1_c, T2_c, f0_c, B1);
    VectorXcd r = m_ab + m_c;
    return r;
}

VectorXcd Three_SSFP_Finite(carrd &flip, const bool spoil,
                            cdbl TR, cdbl Trf, cdbl TE, cdbl ph, cdbl PD,
                            cdbl T1_a, cdbl T2_a,
                            cdbl T1_b, cdbl T2_b,
                            cdbl T1_c, cdbl T2_c,
                            cdbl tau_a, cdbl f_a, cdbl f_c,
                            cdbl f0_a, cdbl f0_b, cdbl f0_c, cdbl B1) {
    double f_ab = 1. - f_c;
    VectorXcd m_ab = Two_SSFP_Finite(flip, spoil, TR, Trf, TE, ph, PD * f_ab, T1_a, T2_a, T1_b, T2_b, tau_a, f_a / f_ab, f0_a, f0_b, B1);
    VectorXcd m_c  = One_SSFP_Finite(flip, spoil, TR, Trf, TE, ph, PD * f_c, T1_c, T2_c, f0_c, B1);
    VectorXcd r = m_ab + m_c;
    return r;
}

} // End namespace QI
