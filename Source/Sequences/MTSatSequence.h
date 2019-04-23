/*
 *  MTSat.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  SPGR / FLASH / FFE Sequences
 *
 */

#ifndef SEQUENCES_MTSAT_H
#define SEQUENCES_MTSAT_H

#include "RFPulse.h"
#include "SequenceBase.h"

namespace QI {

struct MTSatSequence : SequenceBase {
    double         FA, TR, Trf;
    Eigen::ArrayXd sat_f0, sat_angle;
    RFPulse        pulse;
    QI_SEQUENCE_DECLARE(MTSat);
    Eigen::Index size() const override;
};
void from_json(const json &j, MTSatSequence &s);
void to_json(json &j, const MTSatSequence &s);

} // End namespace QI

#endif // SEQUENCES_SPGR_H
