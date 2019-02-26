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

#include <fstream>

#include "JSON.h"
#include "Log.h"
#include "rapidjson/error/en.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/ostreamwrapper.h"
#include "rapidjson/prettywriter.h"

using namespace std::string_literals;

namespace QI {

rapidjson::Document ReadJSON(std::istream &is) {
    rapidjson::IStreamWrapper isw(is);
    rapidjson::Document       doc;
    doc.ParseStream(isw);
    if (doc.HasParseError()) {
        QI::Fail("JSON Error @offset {}: {}",
                 (unsigned)doc.GetErrorOffset(),
                 rapidjson::GetParseError_En(doc.GetParseError()));
    }
    return doc;
}

rapidjson::Document ReadJSON(const std::string &path) {
    std::ifstream ifs(path);
    if (ifs) {
        return ReadJSON(ifs);
    } else {
        QI::Fail("Error reading from JSON file: {}", path);
    }
}

std::ostream &WriteJSON(std::ostream &os, const rapidjson::Document &doc) {
    rapidjson::OStreamWrapper                          osw(os);
    rapidjson::PrettyWriter<rapidjson::OStreamWrapper> writer(osw);
    doc.Accept(writer);
    os << std::endl;
    return os;
}

void WriteJSON(const std::string &path, const rapidjson::Document &doc) {
    std::ofstream ofs(path);
    if (ofs) {
        WriteJSON(ofs, doc);
    } else {
        QI::Fail("Error writing to JSON file: {}", path);
    }
}

rapidjson::Value ArrayToJSON(const Eigen::ArrayXd &              array,
                             rapidjson::Document::AllocatorType &allocator,
                             const double &                      scale) {
    rapidjson::Value json_array(rapidjson::kArrayType);
    for (Eigen::Index i = 0; i < array.rows(); i++) {
        double scaled_value = array[i] * scale;
        json_array.PushBack(scaled_value, allocator);
    }
    return json_array;
}

Eigen::ArrayXd
ArrayFromJSON(const rapidjson::Value &json, const std::string &key, const double &scale) {
    if (!json.HasMember(key))
        QI::Fail("Missing JSON member: {}", key);
    const auto &json_array = json[key];
    if (!json_array.IsArray())
        QI::Fail("JSON member {} is not an array", key);
    Eigen::ArrayXd array(json_array.Size());
    for (rapidjson::SizeType i = 0; i < json_array.Size(); i++) {
        array[i] = json_array[i].GetDouble() * scale;
    }
    return array;
}

const rapidjson::Value &GetMember(const rapidjson::Value &json, const std::string &key) {
    if (!json.HasMember(key))
        QI::Fail("Missing JSON member: {}", key);
    return json[key];
}

} // End namespace QI
