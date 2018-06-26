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

#ifndef SEQUENCES_RFPULSE_H
#define SEQUENCES_RFPULSE_H

#include <string>
#include <iostream>

namespace QI {

struct RFPulse {
    double flip, Trf, intB1, intB1sq;
    std::string name;
};

} // End namespace QI

std::ostream& operator<<(std::ostream &os, const QI::RFPulse &pulse);
std::istream& operator>>(std::istream &is, QI::RFPulse &pulse);

#endif // SEQUENCES_RFPULSE_H
