/*
 *  SequenceGroup.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "SequenceGroup.h"

namespace QI {

size_t SequenceGroup::count() const {
    return sequences.size();
}

const SequenceBase *SequenceGroup::at(const size_t i) const {
    return sequences.at(i);
}

const SequenceBase *SequenceGroup::operator[](const size_t i) const {
    return sequences.at(i);
}

Eigen::Index SequenceGroup::size() const {
    size_t sz = 0;
    for (auto &sig : sequences)
        sz += sig->size();
    return sz;
}

Eigen::ArrayXd SequenceGroup::weights(const double f0) const {
    Eigen::ArrayXd weights(size());
    size_t         start = 0;
    for (auto &sig : sequences) {
        weights.segment(start, sig->size()) = sig->weights(f0);
        start += sig->size();
    }
    return weights;
}

void SequenceGroup::addSequence(const SequenceBase *s) {
    sequences.push_back(s);
}

// SequenceGroup::SequenceGroup(const rapidjson::Value &json) {
//     assert(json.IsArray());
//     for (rapidjson::SizeType i = 0; i < json.Size(); i++) {
//         sequences.push_back(SequenceFromJSON(json[i]));
//     }
// }

// rapidjson::Value SequenceGroup::toJSON(rapidjson::Document::AllocatorType &a) const {
//     rapidjson::Value json(rapidjson::kArrayType);
//     for (const auto &seq : sequences) {
//         auto aval = JSONFromSequence(seq, a);
//         json.PushBack(aval, a);
//     }
//     return json;
// }

} // End namespace QI
