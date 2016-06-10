/*
 *  MPRAGE.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Signals/MPRAGE.h"

using namespace std;
using namespace Eigen;

namespace QI {

VectorXcd One_MPRAGE(cdbl flip, cdbl TR, const int Nseg, const int Nk0, carrd &TI, carrd &TD, cdbl PD, cdbl T1, cdbl B1, cdbl eta) {
    carrd TIs = TI - TR*Nk0; // Adjust TI for k0
    const double M0 = PD;
    const double T1s = 1. / (1./T1 - log(cos(flip * B1))/TR);
    const double M0s = M0 * (1. - exp(-TR/T1)) / (1 - exp(-TR/T1s));

    const double A_1 = M0s*(1 - exp(-(Nseg*TR)/T1s));

    carrd A_2 = M0*(1 - exp(-TD/T1));
    carrd A_3 = M0*(1 - exp(-TIs/T1));
    const double B_1 = exp(-(Nseg*TR)/T1s);
    carrd B_2 = exp(-TD/T1);
    carrd B_3 = -eta*exp(-TIs/T1); // eta is inversion efficency

    carrd A = A_3 + A_2*B_3 + A_1*B_2*B_3;
    carrd B = B_1*B_2*B_3;
    carrd M1 = A / (1. - B);

    VectorXcd M(TI.size());
    M.real() = (M0s + (M1 - M0s)*exp(-(Nk0*TR)/T1s)) * sin(flip * B1);
    M.imag().setZero();
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

} // End namespace QI
