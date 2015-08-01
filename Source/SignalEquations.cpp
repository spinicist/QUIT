/*
 *  SignalEquations.cpp
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright (c) 2011-2013 Tobias Wood.
 *
 *  Based in part on work by Sean Deoni
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "SignalEquations.h"

// There's a bug in matrix log that leads to duplicate symbol definitions if
// this is included in the header file
#include <unsupported/Eigen/MatrixFunctions>

double clamp(double value, double low, double high)
{
	if (value > low) {
		if (value < high) {
			return value;
		} else {
			return high;
		}
	} else {
		return low;
	}
}

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

inline const Matrix3d Relax(const double &T1, const double &T2) {
	Matrix3d R;
	R << 1./T2,     0,     0,
	         0, 1./T2,     0,
	         0,     0, 1./T1;
	return R;
}

inline const Matrix3d InfinitesimalRF(const double &dalpha) {
	Matrix3d A;
	A << 0,      0, -dalpha,
	     0,      0,       0,
	     dalpha, 0,       0;
	return A;
}

inline const Matrix3d OffResonance(const double &inHertz) {
	// Minus signs are this way round to make phase cycling go the right way
	Matrix3d O;
	double dw = inHertz * 2. * M_PI;
	O <<  0, dw, 0,
	    -dw,  0, 0,
		  0,  0, 0;
	return O;
}

inline const Matrix3d Spoiling() {
	// Destroy the x- and y- magnetization
	Matrix3d S = Matrix3d::Zero();
	S(2, 2) = 1.;
	return S;
}

inline const Matrix6d Exchange(const double &k_ab, const double &k_ba) {
	Matrix6d K = Matrix6d::Zero();
	K.block(0,0,3,3).diagonal().setConstant(k_ab);
	K.block(3,3,3,3).diagonal().setConstant(k_ba);
	K.block(3,0,3,3).diagonal().setConstant(-k_ab);
	K.block(0,3,3,3).diagonal().setConstant(-k_ba);
	return K;
}

// Calculate the exchange rates from the residence time and fractions
const void CalcExchange(const double tau_a, const double f_a, double &f_b, double &k_ab, double &k_ba) {
	const double feps = numeric_limits<float>::epsilon(); // Because we read from float files
	f_b = 1.0 - f_a;
	k_ab = 1./tau_a; k_ba = k_ab*f_a/f_b;
	if ((fabs(f_a - 1.) <= feps) || (fabs(f_b - 1.) <= feps)) {
		// Only have 1 component, so no exchange
		k_ab = 0.;
		k_ba = 0.;
	}
}

/******************************************************************************
 * One Component Signals
 *****************************************************************************/
VectorXcd One_MultiEcho(carrd &TE, cdbl PD, cdbl T2) {
	VectorXcd M = VectorXcd::Zero(TE.rows());
	M.real() = PD * (-TE / T2).exp();
	return M;
}

VectorXcd One_SPGR(carrd &flip, cdbl TR, cdbl PD, cdbl T1, cdbl B1) {
	VectorXcd M = VectorXcd::Zero(flip.size());
	ArrayXd sa = (B1 * flip).sin();
	ArrayXd ca = (B1 * flip).cos();
	double expT1 = exp(-TR / T1);
	M.real() = PD * ((1. - expT1) * sa) / (1. - expT1*ca);
	return M;
}

VectorXcd One_SPGR_Echo(carrd &flip, cdbl TR, cdbl TE, cdbl PD, cdbl T1, cdbl T2, cdbl B1) {
    VectorXcd M = VectorXcd::Zero(flip.size());
    ArrayXd sa = (B1 * flip).sin();
    ArrayXd ca = (B1 * flip).cos();
    double e1 = exp(-TR / T1);
    double e2 = exp(-TE / T2); // T2' is incorporated into PD in this formalism
    M.real() = PD * e2 * ((1. - e1) * sa) / (1. - e1*ca);
    return M;
}

VectorXcd One_SSFP(carrd &flip, cdbl TR, cdbl phase,
                   cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1) {
    const Vector3d m0(0., 0., PD);
    const Matrix3d E = (-Relax(T1, T2)*TR).exp();
    const Matrix3d O(AngleAxisd(f0*2.*M_PI*TR, Vector3d::UnitZ()));
    const Matrix3d P(AngleAxisd(phase, Vector3d::UnitZ()));
    MagVector m_e(3, flip.size());
    for (int i = 0; i < flip.size(); i++) {
        const Matrix3d A(AngleAxisd(flip[i] * B1, Vector3d::UnitY()));
        const Vector3d m_minus = (Matrix3d::Identity() - P*O*E*A).partialPivLu().solve((1 - exp(-TR/T1)) * m0);
        m_e.col(i).noalias() = A*m_minus;
    }
    return SigComplex(m_e);
}

VectorXcd One_SSFP_Echo(carrd &flip, cdbl TR, cdbl phase,
                        cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1) {
    double TE = TR / 2;
    const Vector3d m0(0., 0., PD);
    const Matrix3d E = (-Relax(T1, T2)*TR).exp();
    const Matrix3d E_TE = (-Relax(T1, T2)*TE).exp();
    const Matrix3d O(AngleAxisd(f0*2.*M_PI*TR, Vector3d::UnitZ()));
    const Matrix3d O_TE(AngleAxisd(f0*M_PI*TR, Vector3d::UnitZ()));
    const Matrix3d P(AngleAxisd(phase, Vector3d::UnitZ()));
    MagVector m_e(3, flip.size());
    for (int i = 0; i < flip.size(); i++) {
        const Matrix3d A(AngleAxisd(flip[i] * B1, Vector3d::UnitY()));
        const Vector3d m_minus = (Matrix3d::Identity() - P*O*E*A).partialPivLu().solve((1 - exp(-TR/T1)) * m0);
        m_e.col(i).noalias() = O_TE*E_TE*A*m_minus + (1 - exp(-TE/T1)) * m0;
    }
    return SigComplex(m_e);
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

VectorXcd One_SSFP_Ellipse(carrd &flip, cdbl TR, cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1) {
	double E1 = exp(-TR / T1);
    double E2 = exp(-TR / T2);

    double theta = M_PI * f0 * TR;
    ArrayXd M = PD * sqrt(E2) * (1 - E1)*sin(flip * B1) / (1 - E1*E2*E2-(E1-E2*E2)*cos(flip * B1));

    VectorXcd result(flip.size());
    result.real() = M * cos(theta);
    result.imag() = M * sin(theta);

    return result;
}

VectorXcd MP_RAGE(cdbl flip, cdbl TR, const int N, carrd &TI, cdbl TD,
                  cdbl PD, cdbl T1, cdbl B1) {
	const double M0 = PD;
	const double T1s = 1. / (1./T1 - log(cos(flip * B1))/TR);
	const double M0s = M0 * (1. - exp(-TR/T1)) / (1 - exp(-TR/T1s));

	const double A_1 = M0s*(1 - exp(-(N*TR)/T1s));
	const double A_2 = M0*(1 - exp(-TD/T1));
	carrd A_3 = M0*(1 - exp(-TI/T1));
	const double B_1 = exp(-(N*TR)/T1s);
	const double B_2 = exp(-TD/T1);
	carrd B_3 = -exp(-TI/T1);

	carrd A = A_3 + A_2*B_3 + A_1*B_2*B_3;
	carrd B = B_1*B_2*B_3;
	carrd M1 = A / (1. - B);

	VectorXcd M(TI.size()); M.setZero();
	M.real() = M1 * sin(flip * B1);
	return M;
}

VectorXcd One_AFI(cdbl flip, cdbl TR1, cdbl TR2, cdbl PD, cdbl T1, cdbl B1) {
	VectorXcd M = VectorXcd::Zero(2);
	const double E1 = exp(-TR1 / T1);
	const double E2 = exp(-TR2 / T1);
	const double s = sin(B1 * flip);
	const double c = cos(B1 * flip);
	M.real()[0] = PD * s * (1. - E2 + (1. - E1)*E2*c) / (1. - E1*E2*c*c);
	M.real()[1] = PD * s * (1. - E1 + (1. - E2)*E1*c) / (1. - E1*E2*c*c);
	return M;
}

/******************************************************************************
 * Two Component Signals
 *****************************************************************************/
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

VectorXcd Two_SPGR_Echo(carrd &flip, cdbl TR, cdbl TE,
                        cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b, cdbl tau_a, cdbl f_a, cdbl B1) {
    Matrix2d A, eATR;
    Vector2d M0, Mobs;
    MagVector signal(3, flip.size()); signal.setZero();
    double k_ab, k_ba, f_b;
    CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    M0 << f_a, f_b;
    A << -(1./T1_a) - k_ab,              k_ba,
                      k_ab, -(1./T1_b) - k_ba;
    eATR = (TR*A).exp();
    Matrix2d e2 = Matrix2d::Zero();
    e2(0,0) = exp(-TE/T2_a); // T2' absorbed into PD as it effects both components equally
    e2(1,1) = exp(-TE/T2_b);
    const Vector2d RHS = (Matrix2d::Identity() - eATR) * M0;
    for (int i = 0; i < flip.size(); i++) {
        const double a = flip[i] * B1;
        Mobs = e2 * (Matrix2d::Identity() - eATR*cos(a)).partialPivLu().solve(RHS * sin(a));
        signal(1, i) = PD * Mobs.sum();
    }
    return SigComplex(signal);
}

VectorXcd Two_SSFP(carrd &flip, const double TR, const double phase,
                   cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b,
                   cdbl tau_a, cdbl f_a, cdbl f0_a, cdbl f0_b, cdbl B1) {
	MagVector signal(3, flip.size());
	Matrix6d R = Matrix6d::Zero();
	R.block(0,0,3,3) = Relax(T1_a, T2_a);
	R.block(3,3,3,3) = Relax(T1_b, T2_b);
	Matrix6d O = Matrix6d::Zero();
	O.block(0,0,3,3) = OffResonance(f0_a);
	O.block(3,3,3,3) = OffResonance(f0_b);
	double k_ab, k_ba, f_b;
	CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
	Matrix6d K = Exchange(k_ab, k_ba);
    Matrix6d L = (-TR*(R+O+K)).exp();
	Vector6d M0; M0 << 0., 0., PD * f_a, 0., 0., PD * f_b;
	const Vector6d eyemaM0 = (Matrix6d::Identity() - L) * M0;
	Matrix6d A = Matrix6d::Zero();
	for (int i = 0; i < flip.size(); i++) {
		A.block(0, 0, 3, 3) = Matrix3d(AngleAxisd(flip[i] * B1, Vector3d::UnitY()) * AngleAxisd(phase, Vector3d::UnitZ()));
		A.block(3, 3, 3, 3).noalias() = A.block(0, 0, 3, 3);
        Vector6d MTR = (Matrix6d::Identity() - L * A).partialPivLu().solve(eyemaM0);
		signal.col(i) = SumMC(MTR);
	}
	return SigComplex(signal);
}

VectorXcd Two_SSFP_Echo(carrd &flip, const double TR, const double phase,
                        cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b,
                        cdbl tau_a, cdbl f_a, cdbl f0_a, cdbl f0_b, cdbl B1) {
    MagVector signal(3, flip.size());
    Matrix6d R = Matrix6d::Zero();
    R.block(0,0,3,3) = Relax(T1_a, T2_a);
    R.block(3,3,3,3) = Relax(T1_b, T2_b);
    Matrix6d O = Matrix6d::Zero();
    O.block(0,0,3,3) = OffResonance(f0_a);
    O.block(3,3,3,3) = OffResonance(f0_b);
    double k_ab, k_ba, f_b;
    CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    Matrix6d K = Exchange(k_ab, k_ba);
    Matrix6d L2 = ((-TR/2)*(R+O+K)).exp();
    Matrix6d L = L2*L2;
    Vector6d M0; M0 << 0., 0., PD * f_a, 0., 0., PD * f_b;
    const Vector6d eyemaM0 = (Matrix6d::Identity() - L) * M0;
    Matrix6d A = Matrix6d::Zero();
    for (int i = 0; i < flip.size(); i++) {
        A.block(0, 0, 3, 3) = Matrix3d(AngleAxisd(flip[i] * B1, Vector3d::UnitY()) * AngleAxisd(phase, Vector3d::UnitZ()));
        A.block(3, 3, 3, 3).noalias() = A.block(0, 0, 3, 3);
        Vector6d MTR = L2 * (Matrix6d::Identity() - L * A).partialPivLu().solve(eyemaM0);
        signal.col(i) = SumMC(MTR);
    }
    return SigComplex(signal);
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
		l1 = (-(RpOpK+A)*Trf).exp();
		Vector6d m1 = (RpO + A).partialPivLu().solve(Rm0);
		mp.noalias() = Cm2 + (I - l1*C*l2).partialPivLu().solve((I - l1)*(m1 - Cm2));
		me.noalias() = le*(mp - m2) + m2;				
		theory.col(i) = SumMC(me);
	}
	return SigComplex(theory);
}

/******************************************************************************
 * Three Component
 *****************************************************************************/

VectorXcd Three_SPGR(carrd &flip, cdbl TR, cdbl PD,
                     cdbl T1_a, cdbl T1_b, cdbl T1_c, cdbl tau_a, cdbl f_a, cdbl f_c, cdbl B1) {
	double f_ab = 1. - f_c;
	VectorXcd m_ab = Two_SPGR(flip, TR, PD * f_ab, T1_a, T1_b, tau_a, f_a / f_ab, B1);
	VectorXcd m_c  = One_SPGR(flip, TR, PD * f_c, T1_c, B1);
	VectorXcd r = m_ab + m_c;
	return r;
}

VectorXcd Three_SSFP(carrd &flip, cdbl TR, cdbl phase, cdbl PD,
                     cdbl T1_a, cdbl T2_a,
					 cdbl T1_b, cdbl T2_b,
					 cdbl T1_c, cdbl T2_c,
					 cdbl tau_a, cdbl f_a, cdbl f_c,
					 cdbl f0_a, cdbl f0_b, cdbl f0_c, cdbl B1) {
	double f_ab = 1. - f_c;
	VectorXcd m_ab = Two_SSFP(flip, TR, phase, PD * f_ab, T1_a, T2_a, T1_b, T2_b, tau_a, f_a / f_ab, f0_a, f0_b, B1);
	VectorXcd m_c  = One_SSFP(flip, TR, phase, PD * f_c, T1_c, T2_c, f0_c, B1);
	VectorXcd r = m_ab + m_c;
	return r;
}

VectorXcd Three_SSFP_Echo(carrd &flip, cdbl TR, cdbl phase, cdbl PD,
                          cdbl T1_a, cdbl T2_a,
                          cdbl T1_b, cdbl T2_b,
                          cdbl T1_c, cdbl T2_c,
                          cdbl tau_a, cdbl f_a, cdbl f_c,
                          cdbl f0_a, cdbl f0_b, cdbl f0_c, cdbl B1) {
    double f_ab = 1. - f_c;
    VectorXcd m_ab = Two_SSFP_Echo(flip, TR, phase, PD * f_ab, T1_a, T2_a, T1_b, T2_b, tau_a, f_a / f_ab, f0_a, f0_b, B1);
    VectorXcd m_c  = One_SSFP_Echo(flip, TR, phase, PD * f_c, T1_c, T2_c, f0_c, B1);
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

