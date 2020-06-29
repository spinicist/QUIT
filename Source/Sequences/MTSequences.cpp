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

#include "MTSequences.h"

namespace QI {

Eigen::Index ZSpecSequence::size() const {
    return sat_f0.rows();
}

Eigen::Index MTSatSequence::size() const {
    return 1;
}

void from_json(const json &j, ZSpecSequence &s) {
    j.at("TR").get_to(s.TR);
    j.at("Trf").get_to(s.Trf);
    s.FA        = j.at("FA").get<double>() * M_PI / 180.0;
    s.sat_f0    = ArrayFromJSON(j, "sat_f0", 1.);
    s.sat_angle = ArrayFromJSON(j, "sat_angle", M_PI / 180.0);
    j.at("pulse").get_to(s.pulse);
}

void to_json(json &j, const ZSpecSequence &s) {
    j = json{{"TR", s.TR},
             {"FA", s.FA * 180 / M_PI},
             {"pulse", s.pulse},
             {"sat_f0", s.sat_f0},
             {"sat_angle", s.sat_angle * 180 / M_PI}};
}

void from_json(const json &j, MTSatSequence &s) {
    j.at("TR_PDw").get_to(s.TR_pd);
    j.at("TR_T1w").get_to(s.TR_t1);
    j.at("TR_MTw").get_to(s.TR_mt);
    s.al_pd = j.at("FA_PDw").get<double>() * M_PI / 180;
    s.al_t1 = j.at("FA_T1w").get<double>() * M_PI / 180;
    s.al_mt = j.at("FA_MTw").get<double>() * M_PI / 180;
}

void to_json(json &j, const MTSatSequence &s) {
    j = json{
        {"TR_PDw", s.TR_pd},
        {"TR_T1w", s.TR_t1},
        {"TR_MTw", s.TR_mt},
        {"FA_PDw", s.al_pd * 180 / M_PI},
        {"FA_T1w", s.al_t1 * 180 / M_PI},
        {"FA_MTw", s.al_mt * 180 / M_PI},
    };
}

} // End namespace QI