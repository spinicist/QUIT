/*
 *  OnePoolModel.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef MODEL_ONEPOOL_H
#define MODEL_ONEPOOL_H

#include "ModelBase.h"

namespace QI {
namespace Model {

class OnePool : public ModelBase {
    DECLARE_MODEL_INTERFACE()

    virtual Eigen::VectorXd SSFPEchoMagnitude(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;

    virtual Eigen::VectorXcd MultiEcho(cvecd &params, carrd &TE, cdbl TR) const override;
    virtual Eigen::VectorXcd SPGR(cvecd &params, carrd &a, cdbl TR) const override;
    virtual Eigen::VectorXcd SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const override;
    virtual Eigen::VectorXcd SPGRFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl TE) const override;
    virtual Eigen::VectorXcd MPRAGE(cvecd &params, cdbl a, cdbl TR, const int Nseg, const int Nk0, cdbl eta, cdbl TI, cdbl TD) const override;
    virtual Eigen::VectorXcd AFI(cvecd &params, cdbl a, cdbl TR1, cdbl TR2) const override;
    virtual Eigen::VectorXcd SSFP(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual Eigen::VectorXcd SSFPEcho(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual Eigen::VectorXcd SSFPFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, carrd &phi) const override;
    virtual Eigen::VectorXcd SSFP_GS(cvecd &params, carrd &a, cdbl TR) const override;
};

} // End namespace Model
} // End namespace QI

#endif // MODEL_ONEPOOL_H