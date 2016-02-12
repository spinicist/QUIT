#include "Model.h"
#include "Util.h"

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
VectorXcd Model::SSFPEllipse(cvecd &, carrd &, cdbl) const { QI_EXCEPTION("Function not implemented."); }
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

VectorXcd SCD::SSFPEllipse(cvecd &p, carrd &a, cdbl TR) const {
	return scale(One_SSFP_Ellipse(a, TR, p[0], p[1], p[2], p[3], p[4]));
}

/*****************************************************************************/
/* Two Component DESPOT                                                      */
/*****************************************************************************/

string MCD2::Name() const { return "2C"; }
size_t MCD2::nParameters() const { return 9; }
const vector<string> & MCD2::ParameterNames() const {
	static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "tau_m", "f_m", "f0", "B1"};
	return n;
}

ArrayXXd MCD2::Bounds(const FieldStrength f) const {
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

ArrayXd MCD2::Default(const FieldStrength f) const {
    ArrayXd p(nParameters());
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.4,  0.02, 1.0, 0.08, 0.18, 0.1, 0., 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.45, 0.02, 1.2, 0.08, 0.18, 0.1, 0., 1.0; break;
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
    return scale(Two_SPGR_Echo(a, TR, TE, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]));
}

VectorXcd MCD2::SPGRFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl TE) const {
	return scale(Two_SSFP_Finite(a, true, TR, Trf, TE, 0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]));
}

VectorXcd MCD2::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(Two_SSFP(a, phi, TR, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]));
}

VectorXcd MCD2::SSFPEcho(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(Two_SSFP_Echo(a, phi, TR, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[7], p[8]));
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
const vector<string> & MCD2_NoEx::ParameterNames() const {
    static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "f_m", "f0", "B1"};
    return n;
}

ArrayXXd MCD2_NoEx::Bounds(const FieldStrength f) const {
    size_t nP = nParameters();
    ArrayXXd b(nP, 2);
    switch (f) {
        case FieldStrength::Three: b << 1.0, 1.0, 0.3, 0.650, 0.010, 0.030, 0.5, 1.5, 0.05, 0.165, 0.001, 0.35, 0.0, 0.0, 1.0, 1.0; break;
        case FieldStrength::Seven: b << 1.0, 1.0, 0.4, 0.800, 0.010, 0.040, 0.7, 2.5, 0.05, 0.165, 0.001, 0.35, 0.0, 0.0, 1.0, 1.0; break;
        case FieldStrength::User:  b.setZero(); break;
    }
    return b;
}

ArrayXd MCD2_NoEx::Default(const FieldStrength f) const {
    ArrayXd p(nParameters());
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.5, 0.02, 1.0, 0.080, 0.1, 0., 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.7, 0.02, 1.2, 0.075, 0.1, 0., 1.0; break;
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
    return scale(One_SSFP(a, phi, TR, p[0]*p[5], p[1], p[2], p[6], p[7]) +
                 One_SSFP(a, phi, TR, p[0]*(1-p[5]), p[3], p[4], p[6], p[7]));
}

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
    static vector<string> n{"PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "T1_csf", "T2_csf", "tau_m", "f_m", "f_csf", "f0_m", "f0_ie", "B1"};
    return n;
}

ArrayXXd MCD3_f0::Bounds(const FieldStrength f) const {
    size_t nP = nParameters();
    ArrayXXd b(nP, 2);
    switch (f) {
    case FieldStrength::Three: {
        b.col(0) << 1.0, 0.300, 0.010, 0.9, 0.040, 3.5, 1.0, 0.025, 0.001, 0.001, 0.0, 0.0, 1.0;
        b.col(1) << 1.0, 0.800, 0.030, 1.5, 0.150, 5.0, 3.5, 0.600, 0.35,  0.999, 0.0, 0.0, 1.0;
    } break;
    case FieldStrength::Seven: {
        b.col(0) << 1.0, 0.400, 0.010, 0.800, 0.040, 3.0, 0.5, 0.025, 0.001, 0.001, 0.0, 0.0, 1.0;
        b.col(1) << 1.0, 0.650, 0.040, 2.000, 0.140, 4.5, 1.5, 0.600, 0.35,  0.999, 0.0, 0.0, 1.0;
    } break;
    case FieldStrength::User:  b.setZero(); break;
    }
    return b;
}

ArrayXd MCD3_f0::Default(const FieldStrength f) const {
    ArrayXd p(nParameters());
    switch (f) {
    case FieldStrength::Three: p << 1.0, 0.40, 0.015, 1.0, 0.08, 4.0, 3.0, 0.25, 0.1, 0.1, 0.0, 0.0, 1.0; break;
    case FieldStrength::Seven: p << 1.0, 0.45, 0.015, 1.2, 0.08, 4.0, 3.0, 0.25, 0.1, 0.1, 0.0, 0.0, 1.0; break;
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
    return scale(Three_SPGR_Echo(a, TR, TE, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[11], p[12]));
}

VectorXcd MCD3_f0::SPGRFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, cdbl TE) const {
    return scale(Three_SSFP_Finite(a, true, TR, Trf, TE, 0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[11], p[12]));
}

VectorXcd MCD3_f0::SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(Three_SSFP(a, phi, TR, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[11], p[12]));
}

VectorXcd MCD3_f0::SSFPEcho(cvecd &p, carrd &a, cdbl TR, carrd &phi) const {
    return scale(Three_SSFP_Echo(a, phi, TR, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[11], p[12]));
}

VectorXcd MCD3_f0::SSFPFinite(cvecd &p, carrd &a, cdbl TR, cdbl Trf, carrd &phi) const {
   ArrayXcd s(a.rows() * phi.rows());
     ArrayXcd::Index start = 0;
     for (ArrayXcd::Index i = 0; i < phi.rows(); i++) {
         s.segment(start, a.rows()) = Three_SSFP_Finite(a, false, TR, Trf, 0., phi[i], p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[11], p[12]);
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

