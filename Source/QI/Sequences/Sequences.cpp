/*
 *  Sequences.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Sequences/Sequences.h"

namespace QI {

std::shared_ptr<SequenceBase> ReadSequence(std::istream &istr, const bool prompt) {
    if (prompt) std::cout << "Enter sequence type (SPGR/SSFP): " << std::flush;
    std::string type;
    QI::Read(istr, type);
    if (type == "SPGR") {
        return std::make_shared<QI::SPGRSimple>(istr, prompt);
    } else if (type == "SPGR_ECHO") {
        return std::make_shared<QI::SPGREcho>(istr, prompt);
    } else if (type == "SPGR_FINITE") {
        return std::make_shared<QI::SPGRFinite>(istr, prompt);
    } else if (type == "SSFP") {
        return std::make_shared<QI::SSFPSimple>(istr, prompt);
    } else if (type == "SSFP_ECHO") {
        return std::make_shared<QI::SSFPEcho>(istr, prompt);
    } else if (type == "SSFP_ECHO_FLEX") {
        return std::make_shared<QI::SSFPEchoFlex>(istr, prompt);
    } else if (type == "SSFP_FINITE") {
        return std::make_shared<QI::SSFPFinite>(istr, prompt);
    }  else if (type == "SSFP_GS") {
        return std::make_shared<SSFP_GS>(istr, prompt);
    } else if (type == "IRSPGR") {
        return std::make_shared<IRSPGR>(istr, prompt);
    } else if (type == "MPRAGE") {
        return std::make_shared<MPRAGE>(istr, prompt);
    } else if (type == "AFI") {
        return std::make_shared<AFI>(istr, prompt);
    } else if (type == "SPINECHO") {
        return std::make_shared<MultiEcho>(istr, prompt);
    } else {
        QI_EXCEPTION("Unknown sequence type: " << type);
    }
}



} // End namespace QI