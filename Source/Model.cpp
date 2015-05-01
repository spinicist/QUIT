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

string Model::to_string(const Scale &p) {
	static const string sn{"None"}, snm{"Normalised to Mean"};
	switch (p) {
		case Scale::None: return sn;
		case Scale::ToMean: return snm;
	}
}

ArrayXcd Model::scale(const ArrayXcd &s) const {
	ArrayXcd scaled(s.size());
	switch (m_scaling) {
		case Scale::None:   scaled = s; break;
		case Scale::ToMean: scaled = s / s.abs().mean(); break;
	}
	return s;
}

VectorXcd Model::MultiEcho(cvecd &, carrd &) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::SPGR(cvecd &params, carrd &a, cdbl TR) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::SPGRFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl TE) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::MPRAGE(cvecd &params, cdbl a, cdbl TR, const int N, cvecd &TI, cdbl TD) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::SSFP(cvecd &params, carrd &a, cdbl TR, cdbl phi) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::SSFPEllipse(cvecd &params, carrd &a, cdbl TR) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }
VectorXcd Model::SSFPFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl phi) const { throw(logic_error(std::string(__PRETTY_FUNCTION__) + " not implemented.")); }

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

VectorXcd SCD::SSFP(cvecd &p, carrd &a, cdbl TR, cdbl phi) const {
	return scale(One_SSFP(a, TR, phi, p[0], p[1], p[2], p[3], p[4]));
}

VectorXcd SCD::SSFPFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl phi) const {
	return scale(One_SSFP_Finite(a, false, TR, Trf, 0., phi, p[0], p[1], p[2], p[3], p[4]));
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
	static vector<string> n{"PD", "T1_a", "T2_a", "T1_b", "T2_b", "tau_a", "f_a", "f0", "B1"};
	return n;
}

ArrayXXd MCD2::Bounds(const FieldStrength f, cdbl TR) const {
	size_t nP = nParameters();
	ArrayXXd b(nP, 2);
	switch (f) {
		case FieldStrength::Three: b << 1.0, 1.0, 0.1, 0.5, 0.001, 0.030, 0.7, 4.5, 0.05, 0.20, 0.025, 0.60, 0.0, 1.0, -0.5/TR, 0.5/TR, 1.0, 1.0; break;
		case FieldStrength::Seven: b << 1.0, 1.0, 0.1, 0.5, 0.001, 0.025, 1.5, 4.5, 0.04, 0.20, 0.025, 0.60, 0.0, 1.0, -0.5/TR, 0.5/TR, 1.0, 1.0; break;
		case FieldStrength::User:  b.setZero(); break;
	}
	return b;
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

VectorXcd MCD2::SPGRFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl TE) const {
	return scale(Two_SSFP_Finite(a, true, TR, Trf, TE, 0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]));
}

VectorXcd MCD2::SSFP(cvecd &p, carrd &a, cdbl TR, cdbl phi) const {
	return scale(Two_SSFP(a, TR, phi, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]));
}

VectorXcd MCD2::SSFPFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl phi) const {
	return scale(Two_SSFP_Finite(a, false, TR, Trf, 0., phi, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]));
}

/*****************************************************************************/
/* Three Component DESPOT                                                    */
/*****************************************************************************/

string MCD3::Name() const { return "3C"; }
size_t MCD3::nParameters() const { return 12; }
const vector<string> & MCD3::Names() const {
	static vector<string> n{"PD", "T1_a", "T2_a", "T1_b", "T2_b", "T1_c", "T2_c", "tau_a", "f_a", "f_c", "f0", "B1"};
	return n;
}

ArrayXXd MCD3::Bounds(const FieldStrength f, cdbl TR) const {
	size_t nP = nParameters();
	ArrayXXd b(nP, 2);
	switch (f) {
		case FieldStrength::Three: b << 1.0, 1.0, 0.1, 0.5, 0.001, 0.030, 0.7, 2.0, 0.05, 0.20, 3.0, 4.5, 1.5, 2.5, 0.025, 0.60, 0.0, 1.0, 0.0, 1.0, -0.5/TR, 0.5/TR, 1.0, 1.0; break;
		case FieldStrength::Seven: b << 1.0, 1.0, 0.1, 0.5, 0.001, 0.025, 1.5, 2.5, 0.04, 0.20, 3.0, 4.5, 1.5, 2.5, 0.025, 0.60, 0.0, 1.0, 0.0, 1.0, -0.5/TR, 0.5/TR, 1.0, 1.0; break;
		case FieldStrength::User:  b.setZero(); break;
	}
	return b;
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

VectorXcd MCD3::SPGRFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl TE) const {
	return scale(Three_SSFP_Finite(a, true, TR, Trf, TE, 0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]));
}

VectorXcd MCD3::SSFP(cvecd &p, carrd &a, cdbl TR, cdbl phi) const {
	return scale(Three_SSFP(a, TR, phi, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]));
}

VectorXcd MCD3::SSFPFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl phi) const {
	return scale(Three_SSFP_Finite(a, false, TR, Trf, 0., phi, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[10], p[10], p[11]));
}
