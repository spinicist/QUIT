/*
 *  OnePoolModel.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "OnePoolModel.h"
#include "Macro.h"

namespace QI {
namespace Model {

/*****************************************************************************/
/* Single Component DESPOT                                                   */
/*****************************************************************************/

std::string OnePool::Name() const { return "1C"; }
size_t OnePool::nParameters() const { return 5; }
const std::vector<std::string> &OnePool::ParameterNames() const {
    static std::vector<std::string> n{"PD", "T1", "T2", "f0", "B1"};
    return n;
}

Eigen::ArrayXXd OnePool::Bounds(const FieldStrength f) const {
    size_t nP = nParameters();
    Eigen::ArrayXXd b(nP, 2);
    switch (f) {
        case FieldStrength::Three: b << 1.0, 1.0, 0.1, 4.5, 0.010, 2.500, 0.0, 0.0, 1.0, 1.0; break;
        case FieldStrength::Seven: b << 1.0, 1.0, 0.1, 4.5, 0.010, 2.500, 0.0, 0.0, 1.0, 1.0; break;
        case FieldStrength::User:  b.setZero(); break;
    }
    return b;
}

Eigen::ArrayXd OnePool::Default(const FieldStrength /* Unused */) const {
    Eigen::ArrayXd p(5);
    p << 1.0, 1.0, 0.05, 0, 1.0;
    return p;
}

bool OnePool::ValidParameters(cvecd &params) const {
    // Negative T1/T2 makes no sense
    if ((params[1] <= 0.) || (params[2] <= 0.))
        return false;
    else
        return true;
}

Eigen::VectorXcd OnePool::MultiEcho(cvecd &p, carrd &TE, cdbl TR) const {
    return scale(One_MultiEcho(TE, TR, p[0], p[1], p[2]));
}

Eigen::VectorXcd OnePool::SPGR(cvecd &p, carrd &a, cdbl TR) const {
    return scale(One_SPGR(a, TR, p[0], p[1], p[4]));
}

Eigen::VectorXcd OnePool::SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const {
    return scale(One_SPGR_Echo(a, TR, TE, p[0], p[1], p[2], p[3], p[4]));
}
Eigen::VectorXcd OnePool::SPGRFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl TE) const {
    return scale(One_SSFP_Finite(a, true, TR, Trf, TE, 0, p[0], p[1], p[2], p[3], p[4]));
}

Eigen::VectorXcd OnePool::MPRAGE(cvecd &p, cdbl a, cdbl TR, const int Nseg, const int Nk0, cdbl eta, double TI, double TRseg) const {
    return scale(One_MPRAGE(a, TR, Nseg, Nk0, TI, TRseg, p[0], p[1], p[4], eta));
}

Eigen::VectorXcd OnePool::AFI(cvecd &p, cdbl a, cdbl TR1, cdbl TR2) const {
    return scale(One_AFI(a, TR1, TR2, p[0], p[1], p[4]));
}

Eigen::VectorXcd OnePool::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(One_SSFP(a, phi, TR, p[0], p[1], p[2], p[3], p[4]));
}

Eigen::VectorXd OnePool::SSFPEchoMagnitude(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale_mag(One_SSFP_Echo_Magnitude(a, phi, TR, p[0], p[1], p[2], p[3], p[4]));
}

Eigen::VectorXcd OnePool::SSFPEcho(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(One_SSFP_Echo(a, phi, TR, p[0], p[1], p[2], p[3], p[4]));
}

Eigen::VectorXcd OnePool::SSFPFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, carrd &phi) const {
    Eigen::ArrayXcd s(a.rows() * phi.rows());
    Eigen::ArrayXcd::Index start = 0;
    for (Eigen::ArrayXcd::Index i = 0; i < phi.rows(); i++) {
        s.segment(start, a.rows()) = One_SSFP_Finite(a, false, TR, Trf, 0., phi[i], p[0], p[1], p[2], p[3], p[4]);
        start += a.rows();
    }
    return scale(s);
}

Eigen::VectorXcd OnePool::SSFP_GS(cvecd &p, carrd &a, cdbl TR) const {
    return scale(One_SSFP_GS(a, TR, p[0], p[1], p[2], p[3], p[4]));
}

} // End namespace Model
} // End namespace QI