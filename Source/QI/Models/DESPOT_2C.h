/*
 *  DESPOT_2C.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef MODELS_DESPOT_2C_H
#define MODELS_DESPOT_2C_H

#include "QI/Models/Model.h"

namespace QI {

class MCD2 : public Model {
	DECLARE_MODEL_INTERFACE()

    Eigen::VectorXcd SPGR(cvecd &params, carrd &a, cdbl TR) const override;
    Eigen::VectorXcd SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const override;
    Eigen::VectorXcd SPGRFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl TE) const override;
    Eigen::VectorXcd SSFP(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    Eigen::VectorXcd SSFPEcho(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    Eigen::VectorXcd SSFPFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, carrd &phi) const override;
};

class MCD2_NoEx : public Model {
    DECLARE_MODEL_INTERFACE()

    Eigen::VectorXcd SPGR(cvecd &params, carrd &a, cdbl TR) const override;
    Eigen::VectorXcd SSFP(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
};

} // End namespace QI

#endif // MODELS_DESPOT_2C_H