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
VectorXcd One_MultiEcho(carrd &TE, cdbl TR, cdbl PD, cdbl T1, cdbl T2) {
	VectorXcd M = VectorXcd::Zero(TE.rows());
    M.real() = PD * (1 - exp(-TR / T1)) * exp(-TE / T2);
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

VectorXcd One_SSFP_Ellipse(carrd &flip, cdbl TR, cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1) {
    // This is at the echo time
    const double E1 = exp(-TR / T1);
    const double E2 = exp(-TR / T2);

    const double psi_over_2 = M_PI * f0 * TR;
    const ArrayXd alpha = flip * B1;
    const ArrayXcd G = (PD * sqrt(E2) * (1 - E1)*sin(alpha) / (1 - E1*E2*E2-(E1-E2*E2)*cos(alpha))) * polar(1., psi_over_2);

    return G;
}

VectorXcd One_SSFP(carrd &flip, carrd &phi, cdbl TR,
                   cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1) {
    eigen_assert(flip.size() == phi.size())
    const double E1 = exp(-TR / T1);
    const double E2 = exp(-TR / T2);

    const double psi = 2. * M_PI * f0 * TR;
    const ArrayXd alpha = flip * B1;
    const ArrayXd theta = phi + psi;
    // This is not at the echo time
    const ArrayXd d = (1. - E1*E2*E2-(E1-E2*E2)*cos(alpha));
    const ArrayXd G = PD*(1. - E1)*sin(alpha)/d;
    const ArrayXd b = E2*(1. - E1)*(1.+cos(alpha))/d;
    ArrayXcd t(theta.size());
    t.real() = E2 * cos(-theta);
    t.imag() = E2 * sin(theta);
    const ArrayXcd M = G*(1. - t) / (1 - b*cos(theta));
    return M;
}

VectorXcd One_SSFP_Echo(carrd &flip, carrd &phi, cdbl TR,
                        cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1) {
    eigen_assert(flip.size() = phi.size());
    const double E1 = exp(-TR / T1);
    const double E2 = exp(-TR / T2);

    const double psi = 2. * M_PI * f0 * TR;
    const ArrayXd alpha = flip * B1;
    const ArrayXd theta = phi + psi;
    // This is not at the echo time
    const ArrayXd d = (1. - E1*E2*E2-(E1-E2*E2)*cos(alpha));
    const ArrayXcd G = (PD * sqrt(E2) * (1 - E1)*sin(alpha)/d) * polar(1., psi/2.);
    const ArrayXd b = E2*(1. - E1)*(1.+cos(alpha))/d;
    ArrayXcd t(theta.size());
    t.real() = E2 * cos(-theta);
    t.imag() = E2 * sin(theta);
    const ArrayXcd M = G*(1. - t) / (1 - b*cos(theta));
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

VectorXcd One_MPRAGE(cdbl flip, cdbl TR, const int Nseg, const int Nk0, carrd &TI, carrd &TRseg, cdbl PD, cdbl T1, cdbl B1) {
	const double M0 = PD;
	const double T1s = 1. / (1./T1 - log(cos(flip * B1))/TR);
	const double M0s = M0 * (1. - exp(-TR/T1)) / (1 - exp(-TR/T1s));

	const double A_1 = M0s*(1 - exp(-(Nseg*TR)/T1s));

    carrd TD = TRseg - (TI + Nseg*TR);
    carrd A_2 = M0*(1 - exp(-TD/T1));
	carrd A_3 = M0*(1 - exp(-TI/T1));
	const double B_1 = exp(-(Nseg*TR)/T1s);
    carrd B_2 = exp(-TD/T1);
	carrd B_3 = -exp(-TI/T1);

	carrd A = A_3 + A_2*B_3 + A_1*B_2*B_3;
	carrd B = B_1*B_2*B_3;
	carrd M1 = A / (1. - B);

	VectorXcd M(TI.size()); M.setZero();
	M.real() = (M0s + (M1 - M0s)*exp(-Nk0*TR/T1s)) * sin(flip * B1);
	return M;
}

Array2cd One_MP2RAGE(const Array2d &alpha, cdbl TR, const int N, const Array3d &TD,
                  cdbl M0, cdbl T1, cdbl B1, cdbl eta) {
    const double R1 = 1. / T1;
    const Array2d R1s = R1 - log(cos(B1 * alpha))/TR;
    const Array2d M0s = M0 * (1. - exp(-TR*R1)) / (1. - exp(-TR*R1s));
    const double tau = N * TR;

    const Array3d B = exp(-TD*R1);
    const Array3d A = M0*(1. - B);

    const Array2d D = exp(-tau*R1s);
    const Array2d C = M0s*(1. - D);

    Array2d Mm;
    const double denominator = (1 + eta*B[0]*D[0]*B[1]*D[1]*B[2]);
    Mm[0] = (A[0]-eta*B[0]*(A[2]+B[2]*(C[1]+D[1]*(A[1]+B[1]*C[0])))) / denominator;
    Mm[1] = (A[1]+B[1]*(C[0]+D[0]*(A[0]-eta*B[0]*(A[2]+B[2]*C[1])))) / denominator;
    //Mss = -eta*(A[2]+B[2]*(C[1]+D[1]*(A[1]+B[1]*(C[0]+D[0]*A[0])))) / denominator;

    //cout << "denom " << denominator << " Mm " << Mm.transpose() << endl;

    Array2cd Me = Array2cd::Zero();
    Me.real() = Mm * sin(B1 * alpha);
    /*cout << "alpha " << alpha.transpose() << " B1 " << B1 << endl;
    cout << "sin(B1 * alpha) " << sin(B1 * alpha).transpose() << endl;
    cout << "Me " << Me.transpose() << endl;*/
    return Me;
}

Array3cd One_MP3RAGE(const Array3d &alpha, cdbl TR, const int N, const Array4d &TD,
                  cdbl M0, cdbl T1, cdbl B1, cdbl eta) {
    const double R1 = 1. / T1;
    const Array3d R1s = R1 - log(cos(B1 * alpha))/TR;
    const Array3d M0s = M0 * (1. - exp(-TR*R1)) / (1. - exp(-TR*R1s));
    const double tau = N * TR;

    const Array4d B = exp(-TD*R1);
    const Array4d A = M0*(1. - B);

    const Array3d D = exp(-tau*R1s);
    const Array3d C = M0s*(1. - D);

    Array3d Mm;

    const double denominator = (1 + eta*B[0]*D[0]*B[1]*D[1]*B[2]*D[2]*B[3]);
    Mm[0] = (A[0]-eta*B[0]*(A[3]+B[3]*(C[2]+D[2]*(A[2]+B[2]*(C[1]+D[1]*(A[1]+B[1]*C[0])))))) / denominator;
    Mm[1] = (A[1]+B[1]*(C[0]+D[0]*(A[0]-eta*B[0]*(A[3]+B[3]*(C[2]+D[2]*(A[2]+B[2]*C[1])))))) / denominator;
    Mm[2] = (A[2]+B[2]*(C[1]+D[1]*(A[1]+B[1]*(C[0]+D[0]*(A[0]-eta*B[0]*(A[3]+B[3]*C[2])))))) / denominator;

    //cout << "denom " << denominator << " Mm " << Mm.transpose() << endl;

    Array3cd Me = Array3cd::Zero();
    Me.real() = Mm * sin(B1 * alpha);
    /*cout << "alpha " << alpha.transpose() << " B1 " << B1 << endl;
    cout << "sin(B1 * alpha) " << sin(B1 * alpha).transpose() << endl;
    cout << "Me " << Me.transpose() << endl;*/
    return Me;
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

MatrixXd Two_SSFP_Matrix(carrd &flip, carrd &phi, const double TR,
                   cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b,
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
    Vector6d mc;
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
        mc = LHS.partialPivLu().solve(RHS);
        M.col(i) = mc.head(4);
    }
    return M;
}

MatrixXd Two_SSFP_Analytic(carrd &flip, carrd &phi, const double TR,
                   cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b,
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

    const double x0 = E1_a*K4;
    const double x1 = E1_b*K3;
    const double x2 = E1_a*K1;
    carrd ca = cos(alpha);
    carrd sa = sin(alpha);
    carrd sa2 = sa.square();
    carrd cta = cos(theta_a);
    carrd sta = sin(theta_a);
    carrd ctb = cos(theta_b) ;
    carrd stb = sin(theta_b) ;

    carrd x5 = x2 - ca;
    const double x6 = E1_b*K2;
    carrd x7 = x6 + ca;
    carrd x8 = x0*x1 - x5*x7 ;


    carrd x12 = -E2_a*K1*cta;
    carrd x13 = 1/(sa2*x7 + x8*(x12 + ca)) ;
    carrd x14 = x13*x8 ;
    const double x15 = -f_a*x0 + f_b*(-x6 + 1) ;
    carrd x16 = 1/x8 ;
    carrd x17 = x16*sa ;
    carrd x18 = x1*x17 ;

    const double x20 = E2_b*K3 ;
    carrd x21 = sta*x20 ;

    const double x23 = E2_b*K2 ;
    carrd x24 = -ctb*x23 ;

    carrd x26 = x23*stb ;
    const double x27 = E2_a*K4 ;
    carrd x28 = stb*x27 ;
    carrd x29 = E1_b*sa2 + E2_b*cta*x8 ;
    carrd x30 = K3*x13*x29 ;
    carrd x31 = -x26 - x28*x30 ;
    carrd x32 = E1_a*sa2 + E2_a*ctb*x8 ;
    carrd x33 = K4*x13*x32 ;
    carrd x34 = x21*x33 + x26 ;
    carrd x35 = 1/(-K4*x16*x30*x32 + x16*(sa2*x5 + x8*(x24 + ca))) ;
    carrd x36 = x34*x35 ;
    carrd x37 = x13*sta.square()*x8 ;
    carrd x38 = E2_a*K1*x13*sta ;
    carrd x39 = -K3*x29*x38 - x21 ;
    carrd x40 = E2_a*E2_b*K1*K3*x37 - cta*x20 - x36*x39 ;
    const double x41 = E2_a*E2_a ;
    carrd x42 = K4*x32*x38 + x28 ;
    carrd x43 = x35*x42 ;
    carrd x44 = 1/((K1*K1)*x37*x41 + x12 - x39*x43 + 1) ;
    carrd x45 = x44*(K1*K4*x14*sta*stb*x41 - ctb*x27 - x31*x43) ;
    carrd x46 = 1/(E2_a*E2_b*K3*K4*x14*sta*stb + x24 - x31*x36 - x40*x45 + 1) ;
    carrd x47 = E1_b*K3*x13*sa ;
    carrd x48 = -x17*x5 + x18*x33 ;
    carrd x49 = x35*x48 ;
    carrd x50 = E2_a*K1*sta ;
    carrd x51 = -x39*x49 + x47*x50 ;
    carrd x52 = x46*(x28*x47 - x31*x49 - x45*x51) ;
    carrd x53 = x44*(-x40*x52 + x51) ;
    carrd x54 = K3*x16*x29 ;
    carrd x55 = x35*(-x34*x52 - x42*x53 + x48) ;
    const double x56 = f_a*(-x2 + 1) - f_b*x1 ;
    carrd x57 = x17*x7 ;
    carrd x58 = x13*x7*sa ;
    carrd x59 = x0*x17 - x33*x57 ;
    carrd x60 = x35*x59 ;
    carrd x61 = -x39*x60 - x50*x58 ;
    carrd x62 = x46*(-x28*x58 - x31*x60 - x45*x61) ;
    carrd x63 = x44*(-x40*x62 + x61) ;
    carrd x64 = x35*(-x34*x62 - x42*x63 + x59) ;

    MatrixXd M(4, flip.size());
    M.row(0) = x14*x15*(x18 - x21*x52 - x50*x53 + x54*x55) + x14*x56*(-x21*x62 - x50*x63 + x54*x64 - x57);
    M.row(1) = x15*x55 + x56*x64;
    M.row(2) = x15*x53 + x56*x63;
    M.row(3) = x15*x52 + x56*x62;
    return M;
}

VectorXcd Two_SSFP(carrd &flip, carrd &phi, const double TR,
                   cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b,
                   cdbl tau_a, cdbl f_a, cdbl f0_a, cdbl f0_b, cdbl B1) {
    eigen_assert(flip.size() == phi.size())

    MatrixXd M = Two_SSFP_Matrix(flip, phi, TR, PD, T1_a, T2_a, T1_b, T2_b, tau_a, f_a, f0_a, f0_b, B1);
    VectorXcd mc(flip.size());
    mc.real() = M.row(0) + M.row(1);
    mc.imag() = M.row(2) + M.row(3);
    return mc;
}

VectorXcd Two_SSFP_Echo(carrd &flip, carrd &phi, const double TR,
                        cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b,
                        cdbl tau_a, cdbl f_a, cdbl f0_a, cdbl f0_b, cdbl B1) {
    eigen_assert(flip.size() == phi.size())

    MatrixXd M = Two_SSFP_Matrix(flip, phi, TR, PD, T1_a, T2_a, T1_b, T2_b, tau_a, f_a, f0_a, f0_b, B1);

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
    mce.real() = Me.row(0) + Me.row(1);
    mce.imag() = Me.row(2) + Me.row(3);
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

VectorXcd Three_SPGR_Echo(carrd &flip, cdbl TR, cdbl TE,
                          cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b, cdbl T1_c, cdbl T2_c, cdbl tau_a, cdbl f_a, cdbl f_c, cdbl B1) {
    double f_ab = 1. - f_c;
    VectorXcd m_ab = Two_SPGR_Echo(flip, TR, TE, PD * f_ab, T1_a, T2_a, T1_b, T2_b, tau_a, f_a / f_ab, B1);
    VectorXcd m_c  = One_SPGR_Echo(flip, TR, TE, PD * f_c, T1_c, T2_c, B1);
    VectorXcd r = m_ab + m_c;
    return r;
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

    