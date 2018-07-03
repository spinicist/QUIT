/*
 *  JSON.cpp
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "JSON.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/ostreamwrapper.h"
#include "rapidjson/prettywriter.h"

namespace QI {

rapidjson::Document ReadJSON(std::istream &is) {
    rapidjson::IStreamWrapper isw(is);
    rapidjson::Document doc;
    doc.ParseStream(isw);
    return doc;
}

std::ostream &WriteJSON(std::ostream &os, const rapidjson::Document &doc) {
    rapidjson::OStreamWrapper osw(os);
    rapidjson::PrettyWriter<rapidjson::OStreamWrapper> writer(osw);
    doc.Accept(writer);
    return os;
}

void AddArray(rapidjson::Value &val, rapidjson::Document::AllocatorType &a, const std::string &key, const Eigen::ArrayXd &array) {
    rapidjson::Value json_array(rapidjson::kArrayType);
    for (Eigen::Index i = 0; i < array.rows(); i++) {
        json_array.PushBack(array[i], a);
    }
    val.AddMember(rapidjson::Value(key, a), json_array, a);
}

} // End namespace QI
