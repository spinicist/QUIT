/*
 *  MTSat.cpp
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

#include "MTSatSequence.h"
#include "Log.h"

namespace QI {

Eigen::Index MTSatSequence::size() const {
    return sat_f0.rows();
}

void from_json(const json &j, MTSatSequence &s) {
    j.at("TR").get_to(s.TR);
    j.at("Trf").get_to(s.Trf);
    s.FA        = j.at("FA").get<double>() * M_PI / 180.0;
    s.sat_f0    = ArrayFromJSON(j, "sat_f0");
    s.sat_angle = ArrayFromJSON(j, "sat_angle", M_PI / 180.0);
    j.at("pulse").get_to(s.pulse);
}

void to_json(json &j, const MTSatSequence &s) {
    j = json{{"TR", s.TR},
             {"FA", s.FA},
             {"pulse", s.pulse},
             {"sat_f0", s.sat_f0},
             {"sat_angle", s.sat_angle}};
}

} // End namespace QI