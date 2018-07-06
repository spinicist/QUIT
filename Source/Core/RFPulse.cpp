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

RFPulse::RFPulse(const rapidjson::Value &json) {
    Trf = json["Trf"].GetDouble();
    intB1 = json["intB1"].GetDouble();
    intB1sq = json["intB1sq"].GetDouble();
    name = json["name"].GetString();
    FAnom = json["FAnom"].GetDouble() * M_PI / 180;
}

rapidjson::Value RFPulse::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json_val(rapidjson::kObjectType);
    json_val.AddMember("Trf", Trf, a);
    json_val.AddMember("intB1", intB1, a);
    json_val.AddMember("intB1sq", intB1sq, a);
    json_val.AddMember("name", name, a);
    json_val.AddMember("FAnom", FAnom * 180 / M_PI, a);
    return json_val;
}

} // End namespace QI
