/*
 *  RFPulse.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  SPGR sequence with saturation pulse at different offsets
 *
 */

#include "RFPulse.h"
#include "CerealMacro.h"

namespace QI {

void RFPulse::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, Trf);
    QI_CLOAD(ar, intB1);
    QI_CLOAD(ar, intB1sq);
    QI_CLOAD(ar, name);
    QI_CLOAD_DEGREES(ar, FAnom);
}

void RFPulse::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, Trf);
    QI_CSAVE(ar, intB1);
    QI_CSAVE(ar, intB1sq);
    QI_CSAVE(ar, name);
    QI_CSAVE_DEGREES(ar, FAnom);
}

} // End namespace QI
