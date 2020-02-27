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

json ReadJSON(std::string const &path) {
    std::ifstream ifs(path);
    if (ifs) {
        return ReadJSON(ifs);
    } else {
        QI::Fail("Error opening file for reading: {}", path);
    }
}

std::ostream &WriteJSON(std::ostream &os, json const &doc) {
    os << doc.dump(2) << std::endl;
    return os;
}

void WriteJSON(std::string const &path, json const &doc) {
    std::ofstream ofs(path);
    if (ofs) {
        WriteJSON(ofs, doc);
    } else {
        QI::Fail("Could not open file for writing: {}", path);
    }
}

Eigen::ArrayXd ArrayFromJSON(json const &val, std::string const &key, double const &scale) {
    std::vector<double> json_array;
    try {
        json_array = val[key].get<std::vector<double>>();
    } catch (std::exception &e) {
        QI::Fail("Error reading from JSON array {}: {}", key, e.what());
    }
    Eigen::ArrayXd array(json_array.size());
    for (size_t i = 0; i < json_array.size(); i++) {
        array[i] = json_array[i] * scale;
    }
    return array;
}

template <typename T> void GetJSON(json const &j, std::string const &key, T &val) {
    try {
        j.at(key).get_to(val);
    } catch (std::exception &e) {
        QI::Fail("Error reading from JSON value {}: {}", key, e.what());
    }
}

template void GetJSON<double>(json const &j, std::string const &key, double &val);
template void GetJSON<int>(json const &j, std::string const &key, int &val);

} // End namespace QI
