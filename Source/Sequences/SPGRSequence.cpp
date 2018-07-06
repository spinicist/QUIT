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

#include "SPGRSequence.h"
#include "JSON.h"
#include "Macro.h"

namespace QI {

Eigen::Index SPGRBase::size() const {
    return FA.rows();
}

Eigen::ArrayXcd SPGRSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SPGR(p, FA, TR);
}

SPGRSequence::SPGRSequence(const rapidjson::Value &json) {
    if (json.IsNull()) QI_FAIL("Could not read sequence: " << name());
    TR = json["TR"].GetDouble();
    FA = ArrayFromJSON(json["FA"], M_PI / 180);
}

rapidjson::Value SPGRSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("TR", TR, a);
    json.AddMember("FA", ArrayToJSON(FA, a, 180 / M_PI), a);
    return json;
}


/*
 * With echo-time correction
 */

Eigen::ArrayXcd SPGREchoSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SPGREcho(p, FA, TR, TE);
}

SPGREchoSequence::SPGREchoSequence(const rapidjson::Value &json) {
    if (json.IsNull()) QI_FAIL("Could not read sequence: " << name());
    TR = json["TR"].GetDouble();
    TE = json["TE"].GetDouble();
    FA = ArrayFromJSON(json["FA"], M_PI / 180);
}

rapidjson::Value SPGREchoSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("TR", TR, a);
    json.AddMember("TE", TE, a);
    json.AddMember("FA", ArrayToJSON(FA, a, 180 / M_PI), a);
    return json;
}

/*
 * With echo-time and finite-pulse corrections
 */

Eigen::ArrayXcd SPGRFiniteSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    return m->SPGRFinite(p, FA, TR, Trf, TE);
}

SPGRFiniteSequence::SPGRFiniteSequence(const rapidjson::Value &json) {
    if (json.IsNull()) QI_FAIL("Could not read sequence: " << name());
    TR = json["TR"].GetDouble();
    TE = json["TE"].GetDouble();
    Trf = json["Trf"].GetDouble();
    FA = ArrayFromJSON(json["FA"], M_PI / 180);
}

rapidjson::Value SPGRFiniteSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("TR", TR, a);
    json.AddMember("TE", TE, a);
    json.AddMember("Trf", Trf, a);
    json.AddMember("FA", ArrayToJSON(FA, a, 180 / M_PI), a);
    return json;
}

} // End namespace QI