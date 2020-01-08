/*
 *  JSON.h
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_JSON_H
#define QUIT_JSON_H

#include "nlohmann/json.hpp"
#include <Eigen/Core>

using nlohmann::json;

namespace QI {

json          ReadJSON(std::istream &is);
json          ReadJSON(const std::string &path);
std::ostream &WriteJSON(std::ostream &os, const json &doc);
void          WriteJSON(const std::string &path, const json &doc);

Eigen::ArrayXd  ArrayFromJSON(const json &json, const std::string &key, const double &scale = 1);
Eigen::ArrayXcd CArrayFromJSON(const json &json, const std::string &key, const double &scale = 1);

} // End namespace QI

#endif // QUIT_JSON_H