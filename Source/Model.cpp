#include "Model.h"

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
ArrayXcd Model::scale(const ArrayXcd &s) const {
	if (m_scale_to_mean) {
		ArrayXcd scaled = s / s.abs().mean();
		return scaled;
	} else {
		return s;
	}
}

VectorXcd Model::MultiEcho(cvecd &, carrd &) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::SPGR(cvecd &params, carrd &a, cdbl TR) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::SPGRFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl TE) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::MPRAGE(cvecd &params, cdbl a, cdbl TR, const int N, cvecd &TI, cdbl TD) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::AFI(cvecd &p, cdbl a, cdbl TR1, cdbl TR2) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::SSFP(cvecd &params, carrd &a, cdbl TR, carrd &phi) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::SSFPEcho(cvecd &params, carrd &a, cdbl TR, carrd &phi) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::SSFPEllipse(cvecd &params, carrd &a, cdbl TR) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::SSFPFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, carrd &phi) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }

/*****************************************************************************/
/* Single Component DESPOT                                                   */
/*****************************************************************************/

string SCD::Name() const { return "1C"; }
size_t SCD::nParameters() const { return 5; }
const vector<string> &SCD::Names() const {
	static vector<string> n{"PD", "T1", "T2", "f0", "B1"};
	return n;
}

ArrayXXd SCD::Bounds(const FieldStrength f, cdbl TR) const {
	size_t nP = nParameters();
	ArrayXXd b(nP, 2);
	switch (f) {
		case FieldStrength::Three: b << 1.0, 1.0, 0.1, 4.5, 0.010, 2.500, -0.5/TR, 0.5/TR, 1.0, 1.0; break;
		case FieldStrength::Seven: b << 1.0, 1.0, 0.1, 4.5, 0.010, 2.500, -0.5/TR, 0.5/TR, 1.0, 1.0; break;
		case FieldStrength::User:  b.setZero(); break;
	}
	return b;
}

ArrayXd SCD::Start(const FieldStrength f, const double T1, const double T2) const {
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

VectorXcd SCD::MultiEcho(cvecd &p, carrd &TE) const {
	return scale(One_MultiEcho(TE, p[0], p[2]));
}

VectorXcd SCD::SPGR(cvecd &p, carrd &a, cdbl TR) const {
	return scale(One_SPGR(a, TR, p[0], p[1], p[4]));
}

VectorXcd SCD::SPGRFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl TE) const {
	return scale(One_SSFP_Finite(a, true, TR, Trf, TE, 0, p[0], p[1], p[2], p[3], p[4]));
}

VectorXcd SCD::MPRAGE(cvecd &p, cdbl a, cdbl TR, const int N, cvecd &TI, cdbl TD) const {
	return scale(MP_RAGE(a, TR, N, TI, TD, p[0], p[1], p[4]));
}

VectorXcd SCD::AFI(cvecd &p, cdbl a, cdbl TR1, cdbl TR2) const {
	return scale(One_AFI(a, TR1, TR2, p[0], p[1], p[4]));
}

VectorXcd SCD::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    ArrayXcd s(a.rows() * phi.rows());
    ArrayXcd::Index start = 0;
    for (ArrayXcd::Index i = 0; i < phi.rows(); i++) {
        s.segment(start, a.rows()) = One_SSFP(a, TR, phi[i], p[0], p[1], p[2], p[3], p[4]);
        start += a.rows();
    }
    return scale(s);
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

VectorXcd SCD::SSFPEllipse(cvecd &p, carrd &a, cdbl TR) const {
	return scale(One_SSFP_Ellipse(a, TR, p[0], p[1], p[2], p[3], p[4]));
}

/*****************************************************************************/
/* Two Component DESPOT                                                      */
/*****************************************************************************/

string MCD2::Name() const { return "2C"; }
size_t MCD2::nParameters() const { return 9; }
const vector<string> & MCD2::Names() const {
	static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "tau_m", "f_m", "f0", "B1"};
	return n;
}

ArrayXXd MCD2::Bounds(const FieldStrength f, cdbl TR) const {
	size_t nP = nParameters();
	ArrayXXd b(nP, 2);
	switch (f) {
        case FieldStrength::Three: b << 1.0, 1.0, 0.3, 0.65, 0.001, 0.030, 0.5, 1.5, 0.05, 0.165, 0.025, 0.6, 0.001, 0.35, -0.5/TR, 0.5/TR, 1.0, 1.0; break;
        case FieldStrength::Seven: b << 1.0, 1.0, 0.4, 0.8, 0.001, 0.040, 0.7, 2.5, 0.05, 0.165, 0.025, 0.6, 0.001, 0.35, -0.5/TR, 0.5/TR, 1.0, 1.0; break;
		case FieldStrength::User:  b.setZero(); break;
	}
	return b;
}

ArrayXd MCD2::Start(const FieldStrength f, const double T1, const double T2) const {
    ArrayXd p(nParameters());
    const double T2_T1 = T2/T1;
    const double f_m = std::max(0., std::min(0.25, (80. - T2)/240.));
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.5, 0.02, 1.0, 0.08, 0.18, f_m, 0., 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.45, 0.02, 1.1, 0.075, 0.18, f_m, 0., 1.0; break;
    case FieldStrength::User:  p.setZero(); break;
    }
    return p;
}

bool MCD2::ValidParameters(cvecd &params) const {
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

VectorXcd MCD2::SPGR(cvecd &p, carrd &a, cdbl TR) const {
	return scale(Two_SPGR(a, TR, p[0], p[1], p[3], p[5], p[6], p[8]));
}

VectorXcd MCD2::SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const {
    return scale(Two_SPGR_Echo(a, TR, TE, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[8]));
}

VectorXcd MCD2::SPGRFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl TE) const {
	return scale(Two_SSFP_Finite(a, true, TR, Trf, TE, 0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]));
}

VectorXcd MCD2::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    ArrayXcd s(a.rows() * phi.rows());
    ArrayXcd::Index start = 0;
    for (ArrayXcd::Index i = 0; i < phi.rows(); i++) {
        s.segment(start, a.rows()) = Two_SSFP(a, TR, phi[i], p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]);
        start += a.rows();
    }
    return scale(s);
}

VectorXcd MCD2::SSFPEcho(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    ArrayXcd s(a.rows() * phi.rows());
    ArrayXcd::Index start = 0;
    for (ArrayXcd::Index i = 0; i < phi.rows(); i++) {
        s.segment(start, a.rows()) = Two_SSFP_Echo(a, TR, phi[i], p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]);
        start += a.rows();
    }
    return scale(s);
}

VectorXcd MCD2::SSFPFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, carrd &phi) const {
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

string MCD2_NoEx::Name() const { return "2C_NoEx"; }
size_t MCD2_NoEx::nParameters() const { return 8; }
const vector<string> & MCD2_NoEx::Names() const {
    static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "f_m", "f0", "B1"};
    return n;
}

ArrayXXd MCD2_NoEx::Bounds(const FieldStrength f, cdbl TR) const {
    size_t nP = nParameters();
    ArrayXXd b(nP, 2);
    switch (f) {
        case FieldStrength::Three: b << 1.0, 1.0, 0.3, 0.65, 0.010, 0.030, 0.5, 1.5, 0.05, 0.165, 0.001, 0.35, -0.5/TR, 0.5/TR, 1.0, 1.0; break;
        case FieldStrength::Seven: b << 1.0, 1.0, 0.4, 0.8, 0.010, 0.040, 0.7, 2.5, 0.05, 0.165, 0.001, 0.35, -0.5/TR, 0.5/TR, 1.0, 1.0; break;
        case FieldStrength::User:  b.setZero(); break;
    }
    return b;
}

ArrayXd MCD2_NoEx::Start(const FieldStrength f, const double T1, const double T2) const {
    ArrayXd p(nParameters());
    const double T2_T1 = T2/T1;
    const double f_m = std::max(0., std::min(0.25, (80. - T2)/240.));
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.5, 0.02, 1.0, 0.08, f_m, 0., 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.7, 0.02, 1.2, 0.075, f_m, 0., 1.0; break;
    case FieldStrength::User:  p.setZero(); break;
    }
    return p;
}

bool MCD2_NoEx::ValidParameters(cvecd &params) const {
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

VectorXcd MCD2_NoEx::SPGR(cvecd &p, carrd &a, cdbl TR) const {
    return scale(One_SPGR(a, TR, p[0]*p[5], p[1], p[7]) +
                 One_SPGR(a, TR, p[0]*(1-p[5]), p[3], p[7]));
}

VectorXcd MCD2_NoEx::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    ArrayXcd s(a.rows() * phi.rows());
    ArrayXcd::Index start = 0;
    for (ArrayXcd::Index i = 0; i < phi.rows(); i++) {
        s.segment(start, a.rows()) = One_SSFP(a, TR, phi[i], p[0]*p[5], p[1], p[2], p[6], p[7]) +
                                     One_SSFP(a, TR, phi[i], p[0]*(1-p[5]), p[3], p[4], p[6], p[7]);
        start += a.rows();
    }
    return scale(s);
}

/*****************************************************************************/
/* Three Component DESPOT                                                    */
/*****************************************************************************/

string MCD3::Name() const { return "3C"; }
size_t MCD3::nParameters() const { return 12; }
const vector<string> & MCD3::Names() const {
	static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "T1_csf", "T2_csf", "tau_m", "f_m", "f_csf", "f0", "B1"};
	return n;
}

ArrayXXd MCD3::Bounds(const FieldStrength f, cdbl TR) const {
	size_t nP = nParameters();
	ArrayXXd b(nP, 2);
	switch (f) {
    case FieldStrength::Three: {
        b.col(0) << 1.0, 0.300, 0.010, 0.9, 0.040, 3.5, 1.0, 0.025, 0.001, 0.001, -0.5/TR, 1.0;
        b.col(1) << 1.0, 0.800, 0.030, 1.5, 0.150, 5.0, 3.5, 0.600, 0.35,  0.999, 0.5/TR, 1.0;
    } break;
    case FieldStrength::Seven: {
        b.col(0) << 1.0, 0.400, 0.010, 0.800, 0.040, 3.0, 0.5, 0.025, 0.001, 0.001, -0.5/TR, 1.0;
        b.col(1) << 1.0, 0.650, 0.040, 2.000, 0.140, 4.5, 1.5, 0.600, 0.35,  0.999, 0.5/TR, 1.0;
    } break;
    case FieldStrength::User:  b.setZero(); break;
	}
	return b;
}

ArrayXd MCD3::Start(const FieldStrength f, const double T1, const double T2) const {
    ArrayXd p(nParameters());
    const double f_csf = std::max(0.001, std::min(0.5, ((T2/T1) - 0.05) / 0.4));
    const double f_m = 0.15;
    const double T2_ie = std::min(T2 * 1.2, 0.15);
    const double T1_ie = std::min(T1 * 1.2, 1.5);
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.40, 0.015, T1_ie, T2_ie, 4.0, 3.0, 0.25, f_m, f_csf, 0., 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.45, 0.015, T1_ie, T2_ie, 4.0, 3.0, 0.25, f_m, f_csf, 0., 1.0; break;
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
    double f_ab = 1. - p[9];
    VectorXcd m_ab = Two_SPGR_Echo(a, TR, TE, p[0] * f_ab, p[1], p[2], p[3], p[4], p[7], p[8] / f_ab, p[11]);
    VectorXcd m_c  = One_SPGR_Echo(a, TR, TE, p[0] * p[9], p[5], p[6], p[11]);
    VectorXcd r = m_ab + m_c;
    return scale(r);
}

VectorXcd MCD3::SPGRFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl TE) const {
	return scale(Three_SSFP_Finite(a, true, TR, Trf, TE, 0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]));
}

VectorXcd MCD3::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    ArrayXcd s(a.rows() * phi.rows());
    ArrayXcd::Index start = 0;
    for (ArrayXcd::Index i = 0; i < phi.rows(); i++) {
        s.segment(start, a.rows()) = Three_SSFP(a, TR, phi[i], p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]);
        start += a.rows();
    }
    return scale(s);
}

VectorXcd MCD3::SSFPEcho(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    ArrayXcd s(a.rows() * phi.rows());
    ArrayXcd::Index start = 0;
    for (ArrayXcd::Index i = 0; i < phi.rows(); i++) {
        s.segment(start, a.rows()) = Three_SSFP_Echo(a, TR, phi[i], p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]);
        start += a.rows();
    }
    return scale(s);
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
/* Three Component DESPOT w/o/ Exchange                                      */
/*****************************************************************************/

string MCD3_NoEx::Name() const { return "3C_NoEx"; }
size_t MCD3_NoEx::nParameters() const { return 11; }
const vector<string> & MCD3_NoEx::Names() const {
    static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "T1_csf", "T2_csf", "f_m", "f_csf", "f0", "B1"};
    return n;
}

ArrayXXd MCD3_NoEx::Bounds(const FieldStrength f, cdbl TR) const {
    size_t nP = nParameters();
    ArrayXXd b(nP, 2);
    switch (f) {
    case FieldStrength::Three: {
        b.col(0) << 1.0, 0.200, 0.005, 0.700, 0.050, 3.5, 1.0, 0.001, 0.001, -0.5/TR, 1.0;
        b.col(1) << 1.0, 0.700, 0.040, 2.000, 0.200, 5.0, 3.5, 0.500,  0.999, 0.5/TR, 1.0;
    } break;
    case FieldStrength::Seven: {
        b.col(0) << 1.0, 0.400, 0.010, 0.800, 0.040, 3.0, 0.5, 0.025, 0.001, 0.001, -0.5/TR, 1.0;
        b.col(1) << 1.0, 0.650, 0.040, 2.000, 0.140, 4.5, 1.5, 0.600, 0.35,  0.999, 0.5/TR, 1.0;
    } break;
    case FieldStrength::User:  b.setZero(); break;
    }
    return b;
}

ArrayXd MCD3_NoEx::Start(const FieldStrength f, const double T1, const double T2) const {
    ArrayXd p(nParameters());
    const double f_csf = std::max(0.001, std::min(0.5, ((T2/T1) - 0.05) / 0.4));
    const double f_m = 0.15;
    const double T2_ie = std::min(T2 * 1.2, 0.15);
    const double T1_ie = std::min(T1 * 1.2, 1.5);
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.45, 0.015, T1_ie, T2_ie, 4.0, 3.0, f_m, f_csf, 0., 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.45, 0.015, T1_ie, T2_ie, 4.0, 3.0, f_m, f_csf, 0., 1.0; break;
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
    return scale(One_SPGR_Echo(a, TR, TE, p[0]*p[7], p[1], p[2], p[10]) +
                 One_SPGR_Echo(a, TR, TE, p[0]*(1-p[7]-p[8]), p[3], p[4], p[10]) +
                 One_SPGR_Echo(a, TR, TE, p[0]*p[8], p[5], p[6], p[10]));
}


VectorXcd MCD3_NoEx::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    ArrayXcd s(a.rows() * phi.rows());
    ArrayXcd::Index start = 0;
    for (ArrayXcd::Index i = 0; i < phi.rows(); i++) {
        s.segment(start, a.rows()) = One_SSFP(a, TR, phi[i], p[0]*p[7], p[1], p[2], p[9], p[10]) +
                                     One_SSFP(a, TR, phi[i], p[0]*(1-p[7]-p[8]), p[3], p[4], p[9], p[10]) +
                                     One_SSFP(a, TR, phi[i], p[0]*p[8], p[5], p[6], p[9], p[10]);
        start += a.rows();
    }
    return scale(s);
}

VectorXcd MCD3_NoEx::SSFPEcho(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    ArrayXcd s(a.rows() * phi.rows());
    ArrayXcd::Index start = 0;
    for (ArrayXcd::Index i = 0; i < phi.rows(); i++) {
        s.segment(start, a.rows()) = One_SSFP_Echo(a, TR, phi[i], p[0]*p[7], p[1], p[2], p[9], p[10]) +
                                     One_SSFP_Echo(a, TR, phi[i], p[0]*(1-p[7]-p[8]), p[3], p[4], p[9], p[10]) +
                                     One_SSFP_Echo(a, TR, phi[i], p[0]*p[8], p[5], p[6], p[9], p[10]);
        start += a.rows();
    }
    return scale(s);
}

