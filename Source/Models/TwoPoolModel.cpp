/*
 *  TwoPoolModel.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "TwoPoolModel.h"

using namespace std;
using namespace Eigen;

namespace QI {
namespace Model {

/*****************************************************************************/
/* Two Component DESPOT                                                      */
/*****************************************************************************/

string TwoPool::Name() const { return "2C"; }
size_t TwoPool::nParameters() const { return 9; }
const vector<string> & TwoPool::ParameterNames() const {
	static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "tau_m", "f_m", "f0", "B1"};
	return n;
}

ArrayXXd TwoPool::Bounds(const FieldStrength f) const {
	size_t nP = nParameters();
	ArrayXXd b(nP, 2);
	switch (f) {
    case FieldStrength::Three: {
        b.col(0) << 1.0, 0.300, 0.010, 0.9, 0.040, 0.025, 0.001, 0.0, 1.0;
        b.col(1) << 1.0, 0.800, 0.030, 1.5, 0.150, 0.600, 0.35,  0.0, 1.0;
    } break;
    case FieldStrength::Seven: {
        b.col(0) << 1.0, 0.400, 0.010, 0.800, 0.040, 0.025, 0.001, 0.0, 1.0;
        b.col(1) << 1.0, 0.650, 0.040, 2.000, 0.140, 0.600, 0.35,  0.0, 1.0;
    } break;
    case FieldStrength::User:  b.setZero(); break;
	}
	return b;
}

ArrayXd TwoPool::Default(const FieldStrength f) const {
    ArrayXd p(nParameters());
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.4,  0.02, 1.0, 0.08, 0.18, 0.1, 0., 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.45, 0.02, 1.2, 0.08, 0.18, 0.1, 0., 1.0; break;
    case FieldStrength::User:  p.setZero(); break;
    }
    return p;
}

bool TwoPool::ValidParameters(cvecd &params) const {
	// Negative T1/T2 makes no sense
	if ((params[1] <= 0.) || (params[2] <= 0.))
		return false;
	else {
		if ((params[1] < params[3]) && (params[2] < params[4]) && (params[6] <= 1.0))
			return true;
		else
			return false;
	}
}

VectorXcd TwoPool::SPGR(cvecd &p, carrd &a, cdbl TR) const {
	return scale(Two_SPGR(a, TR, p[0], p[1], p[3], p[5], p[6], p[8]));
}

VectorXcd TwoPool::SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const {
    return scale(Two_SPGR_Echo(a, TR, TE, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]));
}

VectorXcd TwoPool::SPGRFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl TE) const {
	return scale(Two_SSFP_Finite(a, true, TR, Trf, TE, 0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]));
}

VectorXcd TwoPool::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(Two_SSFP(a, phi, TR, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]));
}

VectorXcd TwoPool::SSFPEcho(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(Two_SSFP_Echo(a, phi, TR, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]));
}

VectorXcd TwoPool::SSFPFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, carrd &phi) const {
     ArrayXcd s(a.rows() * phi.rows());
     ArrayXcd::Index start = 0;
     for (ArrayXcd::Index i = 0; i < phi.rows(); i++) {
         s.segment(start, a.rows()) = Two_SSFP_Finite(a, false, TR, Trf, 0., phi[i], p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]);
         start += a.rows();
    }
    return scale(s);
}

/*****************************************************************************/
/* Two Component DESPOT w/o exchange                                         */
/*****************************************************************************/

string TwoPool_NoExchange::Name() const { return "2C_NoEx"; }
size_t TwoPool_NoExchange::nParameters() const { return 8; }
const vector<string> & TwoPool_NoExchange::ParameterNames() const {
    static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "f_m", "f0", "B1"};
    return n;
}

ArrayXXd TwoPool_NoExchange::Bounds(const FieldStrength f) const {
    size_t nP = nParameters();
    ArrayXXd b(nP, 2);
    switch (f) {
        case FieldStrength::Three: b << 1.0, 1.0, 0.3, 0.650, 0.010, 0.030, 0.5, 1.5, 0.05, 0.165, 0.001, 0.35, 0.0, 0.0, 1.0, 1.0; break;
        case FieldStrength::Seven: b << 1.0, 1.0, 0.4, 0.800, 0.010, 0.040, 0.7, 2.5, 0.05, 0.165, 0.001, 0.35, 0.0, 0.0, 1.0, 1.0; break;
        case FieldStrength::User:  b.setZero(); break;
    }
    return b;
}

ArrayXd TwoPool_NoExchange::Default(const FieldStrength f) const {
    ArrayXd p(nParameters());
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.5, 0.02, 1.0, 0.080, 0.1, 0., 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.7, 0.02, 1.2, 0.075, 0.1, 0., 1.0; break;
    case FieldStrength::User:  p.setZero(); break;
    }
    return p;
}

bool TwoPool_NoExchange::ValidParameters(cvecd &params) const {
    // Negative T1/T2 makes no sense
    if ((params[1] <= 0.) || (params[2] <= 0.))
        return false;
    else {
        if ((params[1] < params[3]) && (params[2] < params[4]) && (params[6] <= 1.0))
            return true;
        else
            return false;
    }
}

VectorXcd TwoPool_NoExchange::SPGR(cvecd &p, carrd &a, cdbl TR) const {
    return scale(One_SPGR(a, TR, p[0]*p[5], p[1], p[7]) +
                 One_SPGR(a, TR, p[0]*(1-p[5]), p[3], p[7]));
}

VectorXcd TwoPool_NoExchange::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(One_SSFP(a, phi, TR, p[0]*p[5], p[1], p[2], p[6], p[7]) +
                 One_SSFP(a, phi, TR, p[0]*(1-p[5]), p[3], p[4], p[6], p[7]));
}

} // End namespace Model
} // End namespace QI