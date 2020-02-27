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
json          ReadJSON(std::string const &path);
std::ostream &WriteJSON(std::ostream &os, json const &doc);
void          WriteJSON(std::string const &path, json const &doc);

Eigen::ArrayXd ArrayFromJSON(json const &json, std::string const &key, double const &scale = 1);
template <typename T> extern void GetJSON(json const &j, std::string const &key, T &val);

} // End namespace QI

#endif // QUIT_JSON_H