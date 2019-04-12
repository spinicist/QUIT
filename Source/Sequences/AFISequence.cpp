/*
 *  AFI.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "AFISequence.h"
#include "Log.h"

namespace QI {

void from_json(const json &j, AFISequence &s) {
    j.at("TR1").get_to(s.TR1);
    j.at("TR2").get_to(s.TR2);
    j.at("FA").get_to(s.FA);
}

void to_json(json &j, const AFISequence &s) {
    j = nlohmann::json{{"TR1", s.TR1}, {"TR2", s.TR2}, {"FA", s.FA * 180. / M_PI}};
}

} // End namespace QI
