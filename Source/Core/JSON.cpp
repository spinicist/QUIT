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
#include <fstream>
#include <istream>

using namespace std::string_literals;

namespace QI {

json ReadJSON(std::istream &is) {
    return json::parse(is);
}

json ReadJSON(const std::string &path) {
    std::ifstream ifs(path);
    if (ifs) {
        return ReadJSON(ifs);
    } else {
        QI::Fail("Error opening file for reading: {}", path);
    }
}

std::ostream &WriteJSON(std::ostream &os, const json &doc) {
    os << doc.dump(2) << std::endl;
    return os;
}

void WriteJSON(const std::string &path, const json &doc) {
    std::ofstream ofs(path);
    if (ofs) {
        WriteJSON(ofs, doc);
    } else {
        QI::Fail("Could not open file for writing: {}", path);
    }
}

Eigen::ArrayXd ArrayFromJSON(const json &json, const std::string &key, const double &scale) {

    const auto &   json_array = json[key].get<std::vector<double>>();
    Eigen::ArrayXd array(json_array.size());
    for (size_t i = 0; i < json_array.size(); i++) {
        array[i] = json_array[i] * scale;
    }
    return array;
}

} // End namespace QI
