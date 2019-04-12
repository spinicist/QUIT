/*
 *  RFPulse.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  SPGR / FLASH / FFE Sequences
 *
 */

#ifndef QI_RFPULSE_H
#define QI_RFPULSE_H

#include "JSON.h"
#include <string>

namespace QI {

struct RFPulse {
    double Trf, p1, p2, bandwidth;
};
void from_json(const json &j, RFPulse &s);
void to_json(json &j, const RFPulse &s);

} // End namespace QI

#endif // SEQUENCES_RFPULSE_H
