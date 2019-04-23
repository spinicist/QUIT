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

namespace QI {

void from_json(const json &j, RFPulse &s) {
    j.at("p1").get_to(s.p1);
    j.at("p2").get_to(s.p2);
    j.at("bandwidth").get_to(s.bandwidth);
}

void to_json(json &j, const RFPulse &s) {
    j = json{{"p1", s.p1}, {"p2", s.p2}, {"bandwidth", s.bandwidth}};
}

} // End namespace QI
