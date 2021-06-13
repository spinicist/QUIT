/*
 *  SSFP.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

// #define QI_DEBUG_BUILD
#include "Macro.h"

#include "Log.h"
#include "SSFPSequence.h"

namespace QI {

Eigen::Index SSFPBase::size() const {
    return FA.rows();
}

Eigen::ArrayXd SSFPSequence::weights(const double f0) const {
    Eigen::ArrayXd offset  = PhaseInc + M_PI * f0 * TR;
    Eigen::ArrayXd weights = 0.75 * (offset / 2).sin().square();
    QI_DB(f0);
    QI_DBVEC(weights);
    return weights;
}

void from_json(const json &j, SSFPSequence &s) {
    j.at("TR").get_to(s.TR);
    s.FA       = ArrayFromJSON(j, "FA", M_PI / 180.0);
    s.PhaseInc = ArrayFromJSON(j, "PhaseInc", M_PI / 180.0);
    if (s.FA.rows() != s.PhaseInc.rows()) {
        QI::Fail("While reading {} number of phase increments {} did not match the number of "
                 "flip angles {}",
                 s.name(),
                 s.PhaseInc.rows(),
                 s.FA.rows());
    }
}

void to_json(json &j, const SSFPSequence &s) {
    j = json{{"TR", s.TR}, {"FA", s.FA}, {"PhaseInc", s.PhaseInc}};
}

Eigen::ArrayXd SSFPFiniteSequence::weights(const double f0) const {
    Eigen::ArrayXd offset  = PhaseInc + M_PI * f0 * TR;
    Eigen::ArrayXd weights = 0.75 * (offset / 2).sin().square();
    return weights;
}

void from_json(const json &j, SSFPFiniteSequence &s) {
    j.at("TR").get_to(s.TR);
    j.at("Trf").get_to(s.Trf);
    s.FA       = ArrayFromJSON(j, "FA", M_PI / 180.0);
    s.PhaseInc = ArrayFromJSON(j, "PhaseInc", M_PI / 180.0);
    if (s.FA.rows() != s.PhaseInc.rows()) {
        QI::Fail("While reading {} number of phase increments {} did not match the number of "
                 "flip angles {}",
                 s.name(),
                 s.PhaseInc.rows(),
                 s.FA.rows());
    }
}

void to_json(json &j, const SSFPFiniteSequence &s) {
    j = json{{"TR", s.TR}, {"Trf", s.Trf}, {"FA", s.FA}, {"PhaseInc", s.PhaseInc}};
}

Eigen::Index SSFPMTSequence::size() const {
    return FA.rows();
}

void from_json(const json &j, SSFPMTSequence &s) {
    s.TR  = ArrayFromJSON(j, "TR", 1.);
    s.Trf = ArrayFromJSON(j, "Trf", 1.);
    s.FA  = ArrayFromJSON(j, "FA", M_PI / 180);
    j.at("pulse").get_to(s.pulse);
    if ((s.TR.rows() != s.Trf.rows()) || (s.TR.rows() != s.FA.rows())) {
        QI::Fail("Parameters had differing lengths (TR={}, FA={}, Trf={}",
                 s.TR.rows(),
                 s.FA.rows(),
                 s.Trf.rows());
    }
}

void to_json(json &j, const SSFPMTSequence &s) {
    j = json{{"TR", s.TR}, {"Trf", s.Trf}, {"FA", s.FA}, {"pulse", s.pulse}};
}

} // End namespace QI
