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
#include "Log.h"

namespace QI {

void from_json(const json &j, CASLSequence &s) {
    j.at("TR").get_to(s.TR);
    j.at("label_time").get_to(s.label_time);
    s.post_label_delay = ArrayFromJSON(j, "post_label_delay");
}

void to_json(json &j, const CASLSequence &s) {
    j = json{{"TR", s.TR}, {"label_time", s.label_time}, {"post_label_delay", s.post_label_delay}};
}

} // End namespace QI
