/*
 *  DESPOT_3C.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef MODELS_DESPOT_3C_H
#define MODELS_DESPOT_3C_H

#include "QI/Models/Model.h"

namespace QI {
    
 class MCD3 : public Model {
	DECLARE_MODEL_INTERFACE()

    virtual Eigen::VectorXcd SPGR(cvecd &params, carrd &a, cdbl TR) const override;
    virtual Eigen::VectorXcd SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const override;
    virtual Eigen::VectorXcd SPGRFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl TE) const override;
    virtual Eigen::VectorXcd SSFP(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual Eigen::VectorXcd SSFPEcho(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual Eigen::VectorXcd SSFPFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, carrd &phi) const override;
};

class MCD3_f0 : public Model {
    DECLARE_MODEL_INTERFACE()

    virtual Eigen::VectorXcd SPGR(cvecd &params, carrd &a, cdbl TR) const override;
    virtual Eigen::VectorXcd SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const override;
    virtual Eigen::VectorXcd SPGRFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl TE) const override;
    virtual Eigen::VectorXcd SSFP(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual Eigen::VectorXcd SSFPEcho(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual Eigen::VectorXcd SSFPFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, carrd &phi) const override;
};

class MCD3_NoEx : public Model {
    DECLARE_MODEL_INTERFACE()

    virtual Eigen::VectorXcd SPGR(cvecd &p, carrd &a, cdbl TR) const override;
    virtual Eigen::VectorXcd SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const override;
    virtual Eigen::VectorXcd SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const override;
    virtual Eigen::VectorXcd SSFPEcho(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
};

} // End namespace QI

#endif // MODELS_DESPOT_3C_H