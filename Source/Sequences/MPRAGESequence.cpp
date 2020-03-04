/*
QI::Fail("Could not read sequence: {}", name());
 *  MPRAGESequence.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Log.h"
#include "MPRAGESequence.h"

namespace QI {

Eigen::Index MPRAGESequence::size() const {
    return 1;
}

void from_json(const json &j, MPRAGESequence &s) {
    j.at("TR").get_to(s.TR);
    j.at("TI").get_to(s.TI);
    j.at("TD").get_to(s.TD);
    j.at("k0").get_to(s.k0);
    j.at("ETL").get_to(s.ETL);
    j.at("eta").get_to(s.eta);
    s.FA = j.at("FA").get<double>() * M_PI / 180.0;
}

void to_json(json &j, const MPRAGESequence &s) {
    j = {{"TR", s.TR},
         {"TI", s.TI},
         {"TD", s.TD},
         {"eta", s.eta},
         {"FA", s.FA * 180.0 / M_PI},
         {"ETL", s.ETL},
         {"k0", s.k0}};
}

/*
 * MP2RAGE
 */

Eigen::Index MP2RAGESequence::size() const {
    return 2;
}

void from_json(const json &j, MP2RAGESequence &s) {
    j.at("TR").get_to(s.TR);
    auto TI     = ArrayFromJSON(j, "TI", 1.);
    auto TRPrep = j.at("TRPrep").get<double>();
    s.FA        = ArrayFromJSON(j, "FA", M_PI / 180.0);
    j.at("SegLength").get_to(s.SegLength);
    j.at("k0").get_to(s.k0);
    s.TD[0] = TI[0] - s.k0 * s.TR;
    s.TD[1] = TI[1] - (TI[0] + s.SegLength * s.TR);
    s.TD[2] = TRPrep - (TI[1] + (s.SegLength - s.k0) * s.TR);
}

void to_json(json &j, const MP2RAGESequence &s) {
    Eigen::Array2d TI{s.TD[0], s.TD[1] + s.SegLength * s.TR + s.TD[0]};
    double         TRPrep = TI[1] + (s.SegLength - s.k0) * s.TR + s.TD[2];

    j = json{{"TR", s.TR},
             {"TRPrep", TRPrep},
             {"TI", TI},
             {"FA", s.FA},
             {"SegLength", s.SegLength},
             {"k0", s.k0}};
}

} // End namespace QI
