/*
 *  DESPOT_3C.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "DESPOT_3C.h"

using namespace std;
using namespace Eigen;

namespace QI {

/*****************************************************************************/
/* Three Component DESPOT                                                    */
/*****************************************************************************/

string MCD3::Name() const { return "3C"; }
size_t MCD3::nParameters() const { return 12; }
const vector<string> & MCD3::ParameterNames() const {
	static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "T1_csf", "T2_csf", "tau_m", "f_m", "f_csf", "f0", "B1"};
	return n;
}

ArrayXXd MCD3::Bounds(const FieldStrength f) const {
	size_t nP = nParameters();
	ArrayXXd b(nP, 2);
	switch (f) {
    case FieldStrength::Three: {
        b.col(0) << 1.0, 0.300, 0.010, 0.9, 0.040, 3.5, 1.0, 0.025, 0.001, 0.001, 0.0, 1.0;
        b.col(1) << 1.0, 0.800, 0.030, 1.5, 0.150, 5.0, 3.5, 0.600, 0.35,  0.999, 0.0, 1.0;
    } break;
    case FieldStrength::Seven: {
        b.col(0) << 1.0, 0.400, 0.010, 0.800, 0.040, 3.0, 0.5, 0.025, 0.001, 0.001, 0.0, 1.0;
        b.col(1) << 1.0, 0.650, 0.040, 2.000, 0.140, 4.5, 1.5, 0.600, 0.35,  0.999, 0.0, 1.0;
    } break;
    case FieldStrength::User:  b.setZero(); break;
	}
	return b;
}

ArrayXd MCD3::Default(const FieldStrength f) const {
    ArrayXd p(nParameters());
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.40, 0.015, 1.0, 0.08, 4.0, 3.0, 0.25, 0.1, 0.1, 0., 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.45, 0.015, 1.2, 0.08, 4.0, 3.0, 0.25, 0.1, 0.1, 0., 1.0; break;
    case FieldStrength::User:  p.setZero(); break;
    }
    return p;
}

bool MCD3::ValidParameters(cvecd &params) const {
	// Negative T1/T2 makes no sense
	if ((params[1] <= 0.) || (params[2] <= 0.))
		return false;
	else {
		if ((params[1] < params[3]) && (params[2] < params[4]) && (params[3] < params[5]) &&
		    (params[4] < params[6]) && ((params[8] + params[9]) <= 1.0))
			return true;
		else
			return false;
	}
}

VectorXcd MCD3::SPGR(cvecd &p, carrd &a, cdbl TR) const {
	return scale(Three_SPGR(a, TR, p[0], p[1], p[3], p[5], p[7], p[8], p[9], p[11]));
}

VectorXcd MCD3::SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const {
    return scale(Three_SPGR_Echo(a, TR, TE, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]));
}

VectorXcd MCD3::SPGRFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl TE) const {
	return scale(Three_SSFP_Finite(a, true, TR, Trf, TE, 0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]));
}

VectorXcd MCD3::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(Three_SSFP(a, phi, TR, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]));
}

VectorXcd MCD3::SSFPEcho(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(Three_SSFP_Echo(a, phi, TR, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]));
}

VectorXcd MCD3::SSFPFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, carrd &phi) const {
   ArrayXcd s(a.rows() * phi.rows());
     ArrayXcd::Index start = 0;
     for (ArrayXcd::Index i = 0; i < phi.rows(); i++) {
         s.segment(start, a.rows()) = Three_SSFP_Finite(a, false, TR, Trf, 0., phi[i], p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]);
         start += a.rows();
     }
     return scale(s);
}

/*****************************************************************************/
/* Three Component DESPOT with separate myelin off-resonance                 */
/*****************************************************************************/

string MCD3_f0::Name() const { return "3C_f0"; }
size_t MCD3_f0::nParameters() const { return 13; }
const vector<string> & MCD3_f0::ParameterNames() const {
    static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "T1_csf", "T2_csf", "tau_m", "f_m", "f_csf", "f0", "f0_m", "B1"};
    return n;
}

ArrayXXd MCD3_f0::Bounds(const FieldStrength f) const {
    size_t nP = nParameters();
    ArrayXXd b(nP, 2);
    switch (f) {
    case FieldStrength::Three: {
        b.col(0) << 1.0, 0.300, 0.010, 0.9, 0.040, 3.5, 1.0, 0.025, 0.001, 0.001, 0.0, -10.0, 1.0;
        b.col(1) << 1.0, 0.800, 0.030, 1.5, 0.150, 5.0, 3.5, 0.600, 0.35,  0.999, 0.0,  10.0, 1.0;
    } break;
    case FieldStrength::Seven: {
        b.col(0) << 1.0, 0.400, 0.010, 0.800, 0.040, 3.0, 0.5, 0.025, 0.001, 0.001, 0.0, -10.0, 1.0;
        b.col(1) << 1.0, 0.650, 0.040, 2.000, 0.140, 4.5, 1.5, 0.600, 0.35,  0.999, 0.0,  10.0, 1.0;
    } break;
    case FieldStrength::User:  b.setZero(); break;
    }
    return b;
}

ArrayXd MCD3_f0::Default(const FieldStrength f) const {
    ArrayXd p(nParameters());
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.40, 0.015, 1.0, 0.08, 4.0, 3.0, 0.25, 0.1, 0.1, 0.0, -10.0, 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.45, 0.015, 1.2, 0.08, 4.0, 3.0, 0.25, 0.1, 0.1, 0.0, -10.0, 1.0; break;
    case FieldStrength::User:  p.setZero(); break;
    }
    return p;
}

bool MCD3_f0::ValidParameters(cvecd &params) const {
    // Negative T1/T2 makes no sense
    if ((params[1] <= 0.) || (params[2] <= 0.))
        return false;
    else {
        if ((params[1] < params[3]) && (params[2] < params[4]) && (params[3] < params[5]) &&
            (params[4] < params[6]) && ((params[8] + params[9]) <= 1.0))
            return true;
        else
            return false;
    }
}

VectorXcd MCD3_f0::SPGR(cvecd &p, carrd &a, cdbl TR) const {
    return scale(Three_SPGR(a, TR, p[0], p[1], p[3], p[5], p[7], p[8], p[9], p[12]));
}

VectorXcd MCD3_f0::SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const {
    return scale(Three_SPGR_Echo(a, TR, TE, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10]+p[11], p[10], p[10], p[12]));
}

VectorXcd MCD3_f0::SPGRFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl TE) const {
    return scale(Three_SSFP_Finite(a, true, TR, Trf, TE, 0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10]+p[11], p[10], p[10], p[12]));
}

VectorXcd MCD3_f0::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(Three_SSFP(a, phi, TR, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10]+p[11], p[10], p[10], p[12]));
}

VectorXcd MCD3_f0::SSFPEcho(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(Three_SSFP_Echo(a, phi, TR, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10]+p[11], p[10], p[10], p[12]));
}

VectorXcd MCD3_f0::SSFPFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, carrd &phi) const {
   ArrayXcd s(a.rows() * phi.rows());
     ArrayXcd::Index start = 0;
     for (ArrayXcd::Index i = 0; i < phi.rows(); i++) {
         s.segment(start, a.rows()) = Three_SSFP_Finite(a, false, TR, Trf, 0., phi[i], p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10]+p[11], p[10], p[10], p[12]);
         start += a.rows();
     }
     return scale(s);
}

/*****************************************************************************/
/* Three Component DESPOT w/o/ Exchange                                      */
/*****************************************************************************/

string MCD3_NoEx::Name() const { return "3C_NoEx"; }
size_t MCD3_NoEx::nParameters() const { return 11; }
const vector<string> & MCD3_NoEx::ParameterNames() const {
    static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "T1_csf", "T2_csf", "f_m", "f_csf", "f0", "B1"};
    return n;
}

ArrayXXd MCD3_NoEx::Bounds(const FieldStrength f) const {
    size_t nP = nParameters();
    ArrayXXd b(nP, 2);
    switch (f) {
    case FieldStrength::Three: {
        b.col(0) << 1.0, 0.200, 0.005, 0.700, 0.050, 3.5, 1.0, 0.001, 0.001, 0.0, 1.0;
        b.col(1) << 1.0, 0.700, 0.040, 2.000, 0.200, 5.0, 3.5, 0.500, 0.999, 0.0, 1.0;
    } break;
    case FieldStrength::Seven: {
        b.col(0) << 1.0, 0.400, 0.010, 0.800, 0.040, 3.0, 0.5, 0.025, 0.001, 0.001, 0.0, 1.0;
        b.col(1) << 1.0, 0.650, 0.040, 2.000, 0.140, 4.5, 1.5, 0.600, 0.35,  0.999, 0.0, 1.0;
    } break;
    case FieldStrength::User:  b.setZero(); break;
    }
    return b;
}

ArrayXd MCD3_NoEx::Default(const FieldStrength f) const {
    ArrayXd p(nParameters());
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.45, 0.015, 1.0, 0.08, 4.0, 3.0, 0.1, 0.1, 0., 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.45, 0.015, 1.2, 0.08, 4.0, 3.0, 0.1, 0.1, 0., 1.0; break;
    case FieldStrength::User:  p.setZero(); break;
    }
    return p;
}

bool MCD3_NoEx::ValidParameters(cvecd &params) const {
    // Negative T1/T2 makes no sense
    if ((params[1] <= 0.) || (params[2] <= 0.))
        return false;
    else {
        if ((params[1] < params[3]) && (params[2] < params[4]) && (params[3] < params[5]) &&
            (params[4] < params[6]) && ((params[8] + params[9]) <= 1.0))
            return true;
        else
            return false;
    }
}

VectorXcd MCD3_NoEx::SPGR(cvecd &p, carrd &a, cdbl TR) const {
    return scale(One_SPGR(a, TR, p[0]*p[7], p[1], p[10]) +
                 One_SPGR(a, TR, p[0]*(1-p[7]-p[8]), p[3], p[10]) +
                 One_SPGR(a, TR, p[0]*p[8], p[5], p[10]));
}

VectorXcd MCD3_NoEx::SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const {
    return scale(One_SPGR_Echo(a, TR, TE, p[0]*p[7], p[1], p[2], p[9], p[10]) +
                 One_SPGR_Echo(a, TR, TE, p[0]*(1-p[7]-p[8]), p[3], p[4], p[9], p[10]) +
                 One_SPGR_Echo(a, TR, TE, p[0]*p[8], p[5], p[6], p[9], p[10]));
}


VectorXcd MCD3_NoEx::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(One_SSFP(a, phi, TR, p[0]*p[7], p[1], p[2], p[9], p[10]) +
                 One_SSFP(a, phi, TR, p[0]*(1-p[7]-p[8]), p[3], p[4], p[9], p[10]) +
                 One_SSFP(a, phi, TR, p[0]*p[8], p[5], p[6], p[9], p[10]));
}

VectorXcd MCD3_NoEx::SSFPEcho(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(One_SSFP_Echo(a, phi, TR, p[0]*p[7], p[1], p[2], p[9], p[10]) +
                 One_SSFP_Echo(a, phi, TR, p[0]*(1-p[7]-p[8]), p[3], p[4], p[9], p[10]) +
                 One_SSFP_Echo(a, phi, TR, p[0]*p[8], p[5], p[6], p[9], p[10]));
}

} // End namespace QI
