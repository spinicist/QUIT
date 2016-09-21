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
    const ArrayXcd G = polar(PD*sqrt(E2), psi/2)*(1 - E1)*sin(alpha)/d;
    return G;
}

VectorXcd One_SSFP(carrd &flip, carrd &phi, cdbl TR,
                   cdbl PD, cdbl T1, cdbl T2, cdbl f0, cdbl B1) {
    eigen_assert(flip.size() == phi.size());
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
    eigen_assert(flip.size() == phi.size());
    const double E1 = exp(-TR / T1);
    const double E2 = exp(-TR / T2);

    const double   psi = 2. * M_PI * f0 * TR;
    const ArrayXd  alpha = flip * B1;
    const ArrayXd  theta = psi + phi;
    const ArrayXd  d = (1. - E1*E2*E2-(E1-E2*E2)*cos(alpha));
    const ArrayXcd G = polar(PD*sqrt(E2), psi/2.)*(1 - E1)*sin(alpha)/d;
    const ArrayXd  b = E2*(1. - E1)*(1.+cos(alpha))/d;
    ArrayXcd et(theta.size());
    et.real() = cos(-theta);
    et.imag() = sin(-theta);
    const ArrayXcd M = G*(1. - E2*et) / (1 - b*cos(theta));
    /*std::cout << "f0: " << f0 << " TR: " << TR << " psi: " << psi << std::endl;
    std::cout << "phi:   " << phi.transpose() << std::endl;
    std::cout << "theta: " << theta.transpose() << std::endl;
    std::cout << "phase: " << arg(M).transpose() << std::endl;*/
    return M;
}

VectorXd One_SSFP_Echo_Magnitude(carrd &flip, carrd &phi, cdbl TR, cdbl M0, cdbl T1, cdbl T2, cdbl f0, cdbl B1) {
    eigen_assert(flip.size() == phi.size());
    const double E1 = exp(-TR / T1);
    const double E2 = exp(-TR / T2);

    const double  psi = 2. * M_PI * f0 * TR;
    const ArrayXd al = flip * B1;
    const ArrayXd th = phi + psi;
    const ArrayXd d = (1. - E1*cos(al))*(1. - E2*cos(th)) - E2*(E1-cos(al))*(E2-cos(th));
    const ArrayXd rtn = (E2*(E2*E2-2*E2*cos(th)+1.)).sqrt();
    const ArrayXd Mxy = M0*(1.-E1)*rtn*sin(al)/d;
    return Mxy;
}

/*
 * For DESPOT2-FM, only includes M0, T2, f0/th derivs for now
 * Note that the f0 deriv is scaled in terms of f0/TR, i.e. there is no multiply by TR
 * This is to match the scaling of that parameter in DESPOT2-FM which improves fitting
 */
MatrixXd One_SSFP_Echo_Derivs(carrd &flip, carrd &phi, cdbl TR, cdbl M0, cdbl T1, cdbl T2, cdbl f0, cdbl B1) {
    eigen_assert(flip.size() == phi.size());
    const double E1 = exp(-TR / T1);
    const double E2 = exp(-TR / T2);
    const double E2sqr = E2*E2;
    const double psi = 2. * M_PI * f0 * TR;
    const ArrayXd al = flip * B1;
    const ArrayXd sa = sin(al);
    const ArrayXd ca = cos(al);
    const ArrayXd th = phi + psi;
    const ArrayXd sth = sin(th);
    const ArrayXd cth = cos(th);
    const ArrayXd d = (1. - E1*ca)*(1. - E2*cth) - E2*(E1-ca)*(E2-cth);
    const ArrayXd dd = (d*d*(E2sqr - 2*E2*cth + 1)); // Denom for dT2, dth
    const ArrayXd n = E2*(E2sqr-2*E2*cth+1.);
    const ArrayXd rtna = n.sqrt()*(1. - E1)*sa;
    
    MatrixXd drv(flip.size(), 3);
    drv.col(0) = rtna/d;
    drv.col(1) = M0*TR*rtna*(2.*n*(E2*(E1 - ca) + (E1 - ca)*(E2 - cth) - (E1*ca - 1.)*cth) - (E2*(E1 - ca)*(E2 - cth) - (E1*ca - 1.)*(E2*cth - 1.))*(3.*E2sqr - 4.*E2*cth + 1.))/(2.*T2*T2*dd);
    drv.col(2) = 2.*M_PI*E2*M0*rtna*sth*(d + (E2sqr - 2.*E2*cth + 1.)*(E1*ca + E1 - ca - 1.))/dd;
    return drv;
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

} // End namespace QI
