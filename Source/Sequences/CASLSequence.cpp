/*
 *  CASLSequence.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "CASLSequence.h"
#include "Macro.h"

namespace QI {

CASLSequence::CASLSequence(const rapidjson::Value& json) {
    if (json.IsNull()) QI_FAIL("Could not read sequence: " << name());
    TR = json["TR"].GetDouble();
    label_time = json["label_time"].GetDouble();
    post_label_delay = ArrayFromJSON(json["post_label_delay"]);
}

rapidjson::Value CASLSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("TR", TR, a);
    json.AddMember("label_time", label_time, a);
    json.AddMember("post_label_delay", ArrayToJSON(post_label_delay, a), a);
    return json;
}

Eigen::ArrayXcd CASLSequence::signal(const std::shared_ptr<QI::Model::ModelBase> /* Unused */, const Eigen::VectorXd & /* Unused */) const {
    QI_FAIL("Not Implemented");
}

} // End namespace QI
