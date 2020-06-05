/*
 *  MultiEchoSequence.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Log.h"
#include "MultiEchoSequence.h"

namespace QI {

/*
 * Base
 */
Eigen::Index MultiEchoSequence::size() const {
    return TE.rows();
}

void from_json(const json &j, MultiEchoSequence &s) {
    j.at("TR").get_to(s.TR);
    if (j.find("TE") != j.end()) {
        s.TE = ArrayFromJSON(j, "TE", 1.);
    } else {
        auto TE1 = j.at("TE1").get<double>();
        auto ESP = j.at("ESP").get<double>();
        auto ETL = j.at("ETL").get<double>();
        s.TE     = Eigen::ArrayXd::LinSpaced(ETL, TE1, TE1 + ESP * (ETL - 1));
    }
}

void to_json(json &j, const MultiEchoSequence &s) {
    j = json{{"TR", s.TR}, {"TE", s.TE}};
}

} // namespace QI
