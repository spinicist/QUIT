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
#include "Macro.h"

namespace QI {

Eigen::Index MTSatSequence::size() const { return sat_f0.rows(); }

MTSatSequence::MTSatSequence(const rapidjson::Value &json) : pulse(json["pulse"]) {
    if (json.IsNull())
        QI_FAIL("Could not read sequence: " << name());
    TR        = GetMember(json, "TR").GetDouble();
    FA        = GetMember(json, "FA").GetDouble() * M_PI / 180;
    sat_f0    = ArrayFromJSON(json, "sat_f0");
    sat_angle = ArrayFromJSON(json, "sat_angle", M_PI / 180);
}

rapidjson::Value MTSatSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("TR", TR, a);
    json.AddMember("FA", FA * 180 / M_PI, a);
    json.AddMember("pulse", pulse.toJSON(a), a);
    json.AddMember("sat_f0", ArrayToJSON(sat_f0, a), a);
    json.AddMember("sat_angle", ArrayToJSON(sat_angle, a, 180 / M_PI), a);
    return json;
}

} // End namespace QI