/*
 *  MTSat.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  SPGR sequence with saturation pulse at different offsets
 *
 */

#include "MTSatSequence.h"
#include "CerealMacro.h"

namespace QI {

Eigen::Index MTSatSequence::size() const {
    return sat_f0.rows();
}

Eigen::ArrayXcd MTSatSequence::signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    if (auto rm = std::dynamic_pointer_cast<Model::Ramani>(m)) {
        return rm->MTSat(p, FA, TR, sat_f0, sat_angle, pulse);
    } else {
        QI_FAIL("Dyanmic pointer cast to Ramani model failed");
    }
}

void MTSatSequence::load(cereal::JSONInputArchive &ar) {
    QI_CLOAD(ar, TR);
    QI_CLOAD(ar, FA);
    QI_CLOAD(ar, pulse);
    QI_CLOAD(ar, sat_f0);
    QI_CLOAD_DEGREES(ar, sat_angle);
}

void MTSatSequence::save(cereal::JSONOutputArchive &ar) const {
    QI_CSAVE(ar, TR);
    QI_CSAVE(ar, FA);
    QI_CSAVE(ar, pulse);
    QI_CSAVE(ar, sat_f0);
    QI_CSAVE_DEGREES(ar, sat_angle);
}

} // End namespace QI