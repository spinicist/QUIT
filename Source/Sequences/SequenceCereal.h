/*
 *  SequenceCereal.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_CEREAL_H
#define SEQUENCES_CEREAL_H

#include <cereal/archives/json.hpp>
#include "SequenceBase.h"
#include "SPGRSequence.h"
#include "Macro.h"

namespace QI {

template<typename Sequence>
Sequence ReadSequence(std::istream &is, std::string name, bool verbose) {
    cereal::JSONInputArchive in_archive(is);
    Sequence sequence;
    in_archive(sequence);

    if (verbose) {
        std::cout << "Read " << name << ": " << std::endl;
        cereal::JSONOutputArchive archive(std::cout);
        archive(sequence);
    }
    
    return sequence;
}

} // End namespace QI

namespace cereal {

void save(cereal::JSONOutputArchive &ar, std::shared_ptr<QI::SequenceBase> const &s);
void load(cereal::JSONInputArchive &ar, std::shared_ptr<QI::SequenceBase> &sb);

} // End namespace cereal

#endif // SEQUENCES_BASE_H
