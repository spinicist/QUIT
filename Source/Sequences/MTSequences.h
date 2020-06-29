#pragma once
/*
 *  MTSequences.h
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

#include "RFPulse.h"
#include "SequenceBase.h"

namespace QI {

struct ZSpecSequence {
    double         FA, TR, Trf;
    Eigen::ArrayXd sat_f0, sat_angle;
    RFPulse        pulse;
    Eigen::Index   size() const;
};

void from_json(const json &j, ZSpecSequence &s);
void to_json(json &j, const ZSpecSequence &s);

struct MTSatSequence {
    double       TR_pd, al_pd, TR_t1, al_t1, TR_mt, al_mt;
    Eigen::Index size() const;
};

void from_json(const json &j, MTSatSequence &s);
void to_json(json &j, const MTSatSequence &s);

} // End namespace QI
