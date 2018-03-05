/*
 *  SequenceCereal.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCE_CEREAL_H
#define SEQUENCE_CEREAL_H

#include <cereal/archives/json.hpp>
#include "SequenceBase.h"
#include "Macro.h"

namespace QI {

template<typename TSeq>
extern TSeq ReadSequence(cereal::JSONInputArchive &in_archive,
                         bool verbose,
                         std::ostream &os = std::cout);

template<typename TSeq>
extern TSeq ReadSequence(std::istream &is, bool verbose, std::ostream &os = std::cout);

} // End namespace QI

namespace cereal {

void save(cereal::JSONOutputArchive &ar, std::shared_ptr<QI::SequenceBase> const &s);
void load(cereal::JSONInputArchive &ar, std::shared_ptr<QI::SequenceBase> &sb);

} // End namespace cereal

#endif // SEQUENCE_CEREAL_H
