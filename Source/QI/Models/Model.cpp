/*
 *  Model.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Models/Model.h"
#include "QI/Util.h"

using namespace std;
using namespace Eigen;

const string to_string(const FieldStrength& f) {
	static const string f3{"3T"}, f7{"7T"}, fu{"User"};
	switch (f) {
		case FieldStrength::Three: return f3;
		case FieldStrength::Seven: return f7;
		case FieldStrength::User: return fu;
	}
}

/*****************************************************************************/
/* Base Class                                                                */
/*****************************************************************************/
ptrdiff_t Model::ParameterIndex(const string &p) const {
    auto it = find(ParameterNames().begin(), ParameterNames().end(), p);
    if (it != ParameterNames().end()) {
        return distance(ParameterNames().begin(), it);
    } else {
        QI_EXCEPTION("Parameter " << p << " does not exist in model " << Name());
    }
}

ArrayXcd Model::scale(const ArrayXcd &s) const {
    //std::cout << __PRETTY_FUNCTION__ << std::endl;
	if (m_scale_to_mean) {
        //std::cout << "Scaling to mean value" << std::endl;
		ArrayXcd scaled = s / s.abs().mean();
		return scaled;
	} else {
        //std::cout << "Not scaling" << std::endl;
		return s;
	}
}

VectorXcd Model::MultiEcho(cvecd &, carrd &, cdbl) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd Model::SPGR(cvecd &, carrd &, cdbl) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd Model::SPGREcho(cvecd &, carrd &, cdbl, cdbl) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd Model::SPGRFinite(cvecd &, carrd &, cdbl, cdbl, cdbl) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd Model::MPRAGE(cvecd &, cdbl, cdbl, const int, const int, cvecd &, carrd &) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd Model::AFI(cvecd &, cdbl, cdbl, cdbl) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd Model::SSFP(cvecd &, carrd &a, cdbl, carrd &) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd Model::SSFPEcho(cvecd &, carrd &, cdbl, carrd &) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd Model::SSFP_GS(cvecd &, carrd &, cdbl) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd Model::SSFPFinite(cvecd &, carrd &, cdbl, cdbl, carrd &) const { QI_EXCEPTION("Function not implemented."); }

/*****************************************************************************/
/* Single Component DESPOT                                                   */
/*****************************************************************************/

string SCD::Name() const { return "1C"; }
size_t SCD::nParameters() const { return 5; }
const vector<string> &SCD::ParameterNames() const {
	static vector<string> n{"PD", "T1", "T2", "f0", "B1"};
	return n;
}

ArrayXXd SCD::Bounds(const FieldStrength f) const {
	size_t nP = nParameters();
	ArrayXXd b(nP, 2);
	switch (f) {
        case FieldStrength::Three: b << 1.0, 1.0, 0.1, 4.5, 0.010, 2.500, 0.0, 0.0, 1.0, 1.0; break;
        case FieldStrength::Seven: b << 1.0, 1.0, 0.1, 4.5, 0.010, 2.500, 0.0, 0.0, 1.0, 1.0; break;
		case FieldStrength::User:  b.setZero(); break;
	}
	return b;
}

ArrayXd SCD::Default(const FieldStrength f) const {
    ArrayXd p(5);
    p << 1.0, 1.0, 0.05, 0, 1.0;
    return p;
}

bool SCD::ValidParameters(cvecd &params) const {
	// Negative T1/T2 makes no sense
	if ((params[1] <= 0.) || (params[2] <= 0.))
		return false;
	else
		return true;
}

VectorXcd SCD::MultiEcho(cvecd &p, carrd &TE, cdbl TR) const {
    return scale(One_MultiEcho(TE, TR, p[0], p[1], p[2]));
}

VectorXcd SCD::SPGR(cvecd &p, carrd &a, cdbl TR) const {
	return scale(One_SPGR(a, TR, p[0], p[1], p[4]));
}

VectorXcd SCD::SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const {
    return scale(One_SPGR_Echo(a, TR, TE, p[0], p[1], p[2], p[3], p[4]));
}
VectorXcd SCD::SPGRFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl TE) const {
	return scale(One_SSFP_Finite(a, true, TR, Trf, TE, 0, p[0], p[1], p[2], p[3], p[4]));
}

VectorXcd SCD::MPRAGE(cvecd &p, cdbl a, cdbl TR, const int Nseg, const int Nk0, cvecd &TI, carrd &TRseg) const {
    return scale(One_MPRAGE(a, TR, Nseg, Nk0, TI, TRseg, p[0], p[1], p[4]));
}

VectorXcd SCD::AFI(cvecd &p, cdbl a, cdbl TR1, cdbl TR2) const {
	return scale(One_AFI(a, TR1, TR2, p[0], p[1], p[4]));
}

VectorXcd SCD::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(One_SSFP(a, phi, TR, p[0], p[1], p[2], p[3], p[4]));
}

VectorXcd SCD::SSFPEcho(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(One_SSFP_Echo(a, phi, TR, p[0], p[1], p[2], p[3], p[4]));
}

VectorXcd SCD::SSFPFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, carrd &phi) const {
    ArrayXcd s(a.rows() * phi.rows());
    ArrayXcd::Index start = 0;
    for (ArrayXcd::Index i = 0; i < phi.rows(); i++) {
        s.segment(start, a.rows()) = One_SSFP_Finite(a, false, TR, Trf, 0., phi[i], p[0], p[1], p[2], p[3], p[4]);
        start += a.rows();
    }
    return scale(s);
}

VectorXcd SCD::SSFP_GS(cvecd &p, carrd &a, cdbl TR) const {
	return scale(One_SSFP_Ellipse(a, TR, p[0], p[1], p[2], p[3], p[4]));
}
