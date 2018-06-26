/*
 *  qMT.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef MODEL_RAMANI_H
#define MODEL_RAMANI_H

#include "ModelBase.h"
#include "Lineshape.h"

namespace QI {
namespace Model {

class Ramani : public ModelBase {
    DECLARE_MODEL_INTERFACE()
protected:
    TLineshape m_lineshape;
public:
    void setLineshape(TLineshape &l);
    Eigen::VectorXcd SPGR_MT(cvecd &p, carrd &satflip, carrd &satf0, cdbl flip, cdbl TR, cdbl Trf) const override;
};

} // End namespace Model
} // End namespace QI

#endif // MODEL_RAMANI_H
