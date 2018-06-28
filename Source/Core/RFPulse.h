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
#include <iostream>
#include <cereal/archives/json.hpp>

namespace QI {

struct RFPulse {
    double FAnom, Trf, intB1, intB1sq;
    std::string name;
    void load(cereal::JSONInputArchive &ar);
    void save(cereal::JSONOutputArchive &ar) const;
};

} // End namespace QI

#endif // SEQUENCES_RFPULSE_H
