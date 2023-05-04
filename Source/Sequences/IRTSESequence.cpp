/*
 *  IRTSESequence.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  Implements signal equation from 
 *  Padormo, F. et al. In vivo T1 mapping of neonatal brain tissue at 64 mT. 
 *  Magnetic Resonance in Med 89, 1016–1025 (2023).
 * 
 * 
 */

#include "Log.h"
#include "IRTSESequence.h"

namespace QI {

/*
 * Base
 */
Eigen::Index IRTSESequence::size() const {
    return TI.rows();
}

void from_json(const json &j, IRTSESequence &s) {
    j.at("ESP").get_to(s.ESP);
    j.at("ETL").get_to(s.ETL); // Number of refocusing pulses
    j.at("TD1").get_to(s.TD1); // Delay after inversion pulse to navigator pulse
    j.at("theta").get_to(s.theta);
    s.TR = ArrayFromJSON(j, "TR", 1.);
    s.TI = ArrayFromJSON(j, "TI", 1.);
    s.Q = ArrayFromJSON(j, "Q", 1.);    // Effect of inversion pulse (-1 for 180, 0 for 90, 1 for 0deg)
    s.theta = s.theta * M_PI / 180.;
    s.TD2 = s.TR - s.TI - s.TD1 - s.ESP * s.ETL;
}

void to_json(json &j, const IRTSESequence &s) {
    j = json{{"TR", s.TR}};
}

} // namespace QI
