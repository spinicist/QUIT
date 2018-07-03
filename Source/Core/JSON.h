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

#include <iostream>
#include <Eigen/Dense>
#include "rapidjson/document.h"

namespace QI {

rapidjson::Document ReadJSON(std::istream &is);
std::ostream &WriteJSON(std::ostream &os, const rapidjson::Document &doc);
void AddArray(rapidjson::Value &, rapidjson::Document::AllocatorType &, const std::string &, const Eigen::ArrayXd &);
} // End namespace QI

#define QI_JSONIFY( object, doc ) \
    auto value = object.jsonify(doc.GetAllocator()); \
    doc.AddMember( #object, value, doc.GetAllocator());

#endif // QUIT_JSON_H