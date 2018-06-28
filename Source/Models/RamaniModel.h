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
#include "RFPulse.h"

namespace QI {
namespace Model {

class Ramani : public ModelBase {
DECLARE_MODEL_INTERFACE()
protected:
    TLineshape m_lineshape;
public:
    void setLineshape(TLineshape &l);
    Eigen::VectorXd MTSat(cvecd &p, cdbl &FA, cdbl &TR, carrd &satf0, carrd &satflip, const RFPulse &pulse) const;
};

} // End namespace Model
} // End namespace QI

#endif // MODEL_RAMANI_H
