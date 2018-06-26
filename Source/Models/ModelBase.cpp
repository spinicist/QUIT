/*
 *  ModelBase.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "ModelBase.h"
#include "Macro.h"

using namespace std;
using namespace Eigen;

namespace QI {
namespace Model {

/*****************************************************************************/
/* Base Class                                                                */
/*****************************************************************************/
ptrdiff_t ModelBase::ParameterIndex(const string &p) const {
    auto it = find(ParameterNames().begin(), ParameterNames().end(), p);
    if (it != ParameterNames().end()) {
        return distance(ParameterNames().begin(), it);
    } else {
        QI_EXCEPTION("Parameter " << p << " does not exist in model " << Name());
    }
}

ArrayXcd ModelBase::scale(const ArrayXcd &s) const {
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

ArrayXd ModelBase::scale_mag(const ArrayXd &s) const {
    //std::cout << __PRETTY_FUNCTION__ << std::endl;
    if (m_scale_to_mean) {
        //std::cout << "Scaling to mean value" << std::endl;
        ArrayXd scaled = s / s.mean();
        return scaled;
    } else {
        //std::cout << "Not scaling" << std::endl;
        return s;
    }
}

VectorXd ModelBase::SSFPEchoMagnitude(cvecd &, carrd &, cdbl, carrd &) const { QI_EXCEPTION("Function not implemented."); }

VectorXcd ModelBase::MultiEcho(cvecd &, carrd &, cdbl) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd ModelBase::SPGR(cvecd &, carrd &, cdbl) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd ModelBase::SPGREcho(cvecd &, carrd &, cdbl, cdbl) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd ModelBase::SPGRFinite(cvecd &, carrd &, cdbl, cdbl, cdbl) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd ModelBase::SPGR_MT(cvecd &, carrd &, carrd &, cdbl, cdbl, cdbl) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd ModelBase::MPRAGE(cvecd &, cdbl, cdbl, const int, const int, cdbl, cdbl, cdbl) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd ModelBase::AFI(cvecd &, cdbl, cdbl, cdbl) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd ModelBase::SSFP(cvecd &, carrd &, cdbl, carrd &) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd ModelBase::SSFPEcho(cvecd &, carrd &, cdbl, carrd &) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd ModelBase::SSFP_GS(cvecd &, carrd &, cdbl) const { QI_EXCEPTION("Function not implemented."); }
VectorXcd ModelBase::SSFPFinite(cvecd &, carrd &, cdbl, cdbl, carrd &) const { QI_EXCEPTION("Function not implemented."); }

} // End namespace Model
} // End namespace QI
