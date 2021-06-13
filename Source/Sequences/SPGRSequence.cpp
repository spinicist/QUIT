/*
 *  SPGR.cpp
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

#include "JSON.h"
#include "Log.h"
#include "SPGRSequence.h"

namespace QI {

SPGRSequence::SPGRSequence(Eigen::ArrayXd const &FA_, double const &TR_) : FA{FA_}, TR{TR_} {}

Eigen::Index SPGRSequence::size() const {
    return FA.rows();
}

void from_json(const json &j, SPGRSequence &s) {
    j.at("TR").get_to(s.TR);
    s.FA = ArrayFromJSON(j, "FA", M_PI / 180.0);
}

void to_json(json &j, const SPGRSequence &s) {
    j = nlohmann::json{{"TR", s.TR}, {"FA", s.FA * 180.0 / M_PI}};
}

Eigen::Index SPGREchoSequence::size() const {
    return FA.rows();
}

void from_json(const json &j, SPGREchoSequence &s) {
    j.at("TR").get_to(s.TR);
    j.at("TE").get_to(s.TE);
    s.FA = ArrayFromJSON(j, "FA", M_PI / 180.0);
}

void to_json(json &j, const SPGREchoSequence &s) {
    j = nlohmann::json{{"TR", s.TR}, {"TE", s.TE}, {"FA", s.FA * 180.0 / M_PI}};
}

// /*
//  * With echo-time and finite-pulse corrections
//  */
// SPGRFiniteSequence::SPGRFiniteSequence(const rapidjson::Value &json) {
//     if (json.IsNull())
//         QI::Fail("Could not read sequence: {}", name());
//     TR  = GetMember(json, "TR").GetDouble();
//     TE  = GetMember(json, "TE").GetDouble();
//     Trf = GetMember(json, "Trf").GetDouble();
//     FA  = ArrayFromJSON(json, "FA", M_PI / 180);
// }

// rapidjson::Value SPGRFiniteSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
//     rapidjson::Value json(rapidjson::kObjectType);
//     json.AddMember("TR", TR, a);
//     json.AddMember("TE", TE, a);
//     json.AddMember("Trf", Trf, a);
//     json.AddMember("FA", ArrayToJSON(FA, a, 180 / M_PI), a);
//     return json;
// }

} // End namespace QI