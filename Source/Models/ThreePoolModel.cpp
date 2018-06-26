/*
 *  ThreePoolModel.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "ThreePoolModel.h"

using namespace std;
using namespace Eigen;

namespace QI {
namespace Model {

/*****************************************************************************/
/* Three Component DESPOT                                                    */
/*****************************************************************************/

string ThreePool::Name() const { return "3C"; }
size_t ThreePool::nParameters() const { return 12; }
const vector<string> & ThreePool::ParameterNames() const {
	static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "T1_csf", "T2_csf", "tau_m", "f_m", "f_csf", "f0", "B1"};
	return n;
}

ArrayXXd ThreePool::Bounds(const FieldStrength f) const {
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

ArrayXd ThreePool::Default(const FieldStrength f) const {
    ArrayXd p(nParameters());
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.40, 0.015, 1.0, 0.08, 4.0, 3.0, 0.25, 0.1, 0.1, 0., 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.45, 0.015, 1.2, 0.08, 4.0, 3.0, 0.25, 0.1, 0.1, 0., 1.0; break;
    case FieldStrength::User:  p.setZero(); break;
    }
    return p;
}

bool ThreePool::ValidParameters(cvecd &params) const {
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

VectorXcd ThreePool::SPGR(cvecd &p, carrd &a, cdbl TR) const {
	return scale(Three_SPGR(a, TR, p[0], p[1], p[3], p[5], p[7], p[8], p[9], p[11]));
}

VectorXcd ThreePool::SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const {
    return scale(Three_SPGR_Echo(a, TR, TE, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]));
}

VectorXcd ThreePool::SPGRFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl TE) const {
	return scale(Three_SSFP_Finite(a, true, TR, Trf, TE, 0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]));
}

VectorXcd ThreePool::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(Three_SSFP(a, phi, TR, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]));
}

VectorXcd ThreePool::SSFPEcho(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(Three_SSFP_Echo(a, phi, TR, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]));
}

VectorXcd ThreePool::SSFPFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, carrd &phi) const {
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

string ThreePool_f0::Name() const { return "3C_f0"; }
size_t ThreePool_f0::nParameters() const { return 13; }
const vector<string> & ThreePool_f0::ParameterNames() const {
    static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "T1_csf", "T2_csf", "tau_m", "f_m", "f_csf", "f0", "f0_m", "B1"};
    return n;
}

ArrayXXd ThreePool_f0::Bounds(const FieldStrength f) const {
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

ArrayXd ThreePool_f0::Default(const FieldStrength f) const {
    ArrayXd p(nParameters());
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.40, 0.015, 1.0, 0.08, 4.0, 3.0, 0.25, 0.1, 0.1, 0.0, -10.0, 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.45, 0.015, 1.2, 0.08, 4.0, 3.0, 0.25, 0.1, 0.1, 0.0, -10.0, 1.0; break;
    case FieldStrength::User:  p.setZero(); break;
    }
    return p;
}

bool ThreePool_f0::ValidParameters(cvecd &params) const {
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

VectorXcd ThreePool_f0::SPGR(cvecd &p, carrd &a, cdbl TR) const {
    return scale(Three_SPGR(a, TR, p[0], p[1], p[3], p[5], p[7], p[8], p[9], p[12]));
}

VectorXcd ThreePool_f0::SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const {
    return scale(Three_SPGR_Echo(a, TR, TE, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10]+p[11], p[10], p[10], p[12]));
}

VectorXcd ThreePool_f0::SPGRFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl TE) const {
    return scale(Three_SSFP_Finite(a, true, TR, Trf, TE, 0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10]+p[11], p[10], p[10], p[12]));
}

VectorXcd ThreePool_f0::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(Three_SSFP(a, phi, TR, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10]+p[11], p[10], p[10], p[12]));
}

VectorXcd ThreePool_f0::SSFPEcho(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(Three_SSFP_Echo(a, phi, TR, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10]+p[11], p[10], p[10], p[12]));
}

VectorXcd ThreePool_f0::SSFPFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, carrd &phi) const {
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

string ThreePool_NoExchange::Name() const { return "3C_NoEx"; }
size_t ThreePool_NoExchange::nParameters() const { return 11; }
const vector<string> & ThreePool_NoExchange::ParameterNames() const {
    static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "T1_csf", "T2_csf", "f_m", "f_csf", "f0", "B1"};
    return n;
}

ArrayXXd ThreePool_NoExchange::Bounds(const FieldStrength f) const {
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

ArrayXd ThreePool_NoExchange::Default(const FieldStrength f) const {
    ArrayXd p(nParameters());
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.45, 0.015, 1.0, 0.08, 4.0, 3.0, 0.1, 0.1, 0., 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.45, 0.015, 1.2, 0.08, 4.0, 3.0, 0.1, 0.1, 0., 1.0; break;
    case FieldStrength::User:  p.setZero(); break;
    }
    return p;
}

bool ThreePool_NoExchange::ValidParameters(cvecd &params) const {
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

VectorXcd ThreePool_NoExchange::SPGR(cvecd &p, carrd &a, cdbl TR) const {
    return scale(One_SPGR(a, TR, p[0]*p[7], p[1], p[10]) +
                 One_SPGR(a, TR, p[0]*(1-p[7]-p[8]), p[3], p[10]) +
                 One_SPGR(a, TR, p[0]*p[8], p[5], p[10]));
}

VectorXcd ThreePool_NoExchange::SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const {
    return scale(One_SPGR_Echo(a, TR, TE, p[0]*p[7], p[1], p[2], p[9], p[10]) +
                 One_SPGR_Echo(a, TR, TE, p[0]*(1-p[7]-p[8]), p[3], p[4], p[9], p[10]) +
                 One_SPGR_Echo(a, TR, TE, p[0]*p[8], p[5], p[6], p[9], p[10]));
}


VectorXcd ThreePool_NoExchange::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(One_SSFP(a, phi, TR, p[0]*p[7], p[1], p[2], p[9], p[10]) +
                 One_SSFP(a, phi, TR, p[0]*(1-p[7]-p[8]), p[3], p[4], p[9], p[10]) +
                 One_SSFP(a, phi, TR, p[0]*p[8], p[5], p[6], p[9], p[10]));
}

VectorXcd ThreePool_NoExchange::SSFPEcho(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(One_SSFP_Echo(a, phi, TR, p[0]*p[7], p[1], p[2], p[9], p[10]) +
                 One_SSFP_Echo(a, phi, TR, p[0]*(1-p[7]-p[8]), p[3], p[4], p[9], p[10]) +
                 One_SSFP_Echo(a, phi, TR, p[0]*p[8], p[5], p[6], p[9], p[10]));
}

} // End namespace Model 
} // End namespace QI
