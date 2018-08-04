/*
 *  ThreePoolModel.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef MODEL_THREEPOOL_H
#define MODEL_THREEPOOL_H

#include "ModelBase.h"

namespace QI {
namespace Model {

class ThreePool : public ModelBase {
    DECLARE_MODEL_INTERFACE()

    virtual Eigen::VectorXcd SPGR(cvecd &params, carrd &a, cdbl TR) const override;
    virtual Eigen::VectorXcd SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const override;
    virtual Eigen::VectorXcd SPGRFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl TE) const override;
    virtual Eigen::VectorXcd SSFP(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual Eigen::VectorXcd SSFPEcho(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual Eigen::VectorXcd SSFPFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, carrd &phi) const override;
};

class ThreePool_f0 : public ModelBase {
    DECLARE_MODEL_INTERFACE()

    virtual Eigen::VectorXcd SPGR(cvecd &params, carrd &a, cdbl TR) const override;
    virtual Eigen::VectorXcd SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const override;
    virtual Eigen::VectorXcd SPGRFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, cdbl TE) const override;
    virtual Eigen::VectorXcd SSFP(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual Eigen::VectorXcd SSFPEcho(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
    virtual Eigen::VectorXcd SSFPFinite(cvecd &params, carrd &a, cdbl TR, cdbl T_rf, carrd &phi) const override;
};

class ThreePool_NoExchange : public ModelBase {
    DECLARE_MODEL_INTERFACE()

    virtual Eigen::VectorXcd SPGR(cvecd &p, carrd &a, cdbl TR) const override;
    virtual Eigen::VectorXcd SPGREcho(cvecd &p, carrd &a, cdbl TR, cdbl TE) const override;
    virtual Eigen::VectorXcd SSFP(cvecd &p, carrd &a, cdbl TR, carrd &phi) const override;
    virtual Eigen::VectorXcd SSFPEcho(cvecd &params, carrd &a, cdbl TR, carrd &phi) const override;
};

} // End namespace Model
} // End namespace QI

#endif // MODEL_THREEPOOL_H