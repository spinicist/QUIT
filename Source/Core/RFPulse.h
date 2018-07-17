/*
 *  RFPulse.h
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

#ifndef QI_RFPULSE_H
#define QI_RFPULSE_H

#include <string>
#include "JSON.h"

namespace QI {

struct RFPulse {
    double FAnom = 0.0, Trf = 0.0, intB1 = 0.0, intB1sq = 0.0;
    std::string name;

    rapidjson::Value toJSON(rapidjson::Document::AllocatorType &a) const;
    RFPulse(const rapidjson::Value &);
    RFPulse() = default;
};

} // End namespace QI

#endif // SEQUENCES_RFPULSE_H
