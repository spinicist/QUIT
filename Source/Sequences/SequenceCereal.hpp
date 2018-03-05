/*
 *  SequenceCereal.hpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCE_CEREAL_HPP
#define SEQUENCE_CEREAL_HPP

#include <cereal/archives/json.hpp>
#include "SequenceBase.h"
#include "Macro.h"

namespace QI {

template<typename TSeq>
TSeq ReadSequence(cereal::JSONInputArchive &in_archive, bool verbose, std::ostream &os) {
    TSeq sequence;
    try {
        in_archive(cereal::make_nvp(sequence.name(), sequence));
    } catch (cereal::RapidJSONException &e) {
        QI_FAIL("Error parsing input for sequence: " << sequence.name() << "\n" << e.what());
    }

    if (verbose) {
        {   // Archives don't fully flush until destruction
            cereal::JSONOutputArchive archive(std::cout);
            archive(cereal::make_nvp(sequence.name(), sequence));
        }
        std::cout << std::endl;
    }
    return sequence;
}

template<typename TSeq>
TSeq ReadSequence(std::istream &is, bool verbose, std::ostream &os) {
    cereal::JSONInputArchive in_archive(is);
    return ReadSequence<TSeq>(in_archive, verbose, os);
}

} // End namespace QI

#endif // SEQUENCE_CEREAL_HPP