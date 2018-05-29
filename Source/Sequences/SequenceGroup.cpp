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

std::shared_ptr<SequenceBase> &SequenceGroup::operator[](const size_t i) {
    return sequences.at(i);
}

Eigen::Index SequenceGroup::size() const {
    size_t sz = 0;
    for (auto& sig : sequences)
        sz += sig->size();
    return sz;
}

Eigen::ArrayXcd SequenceGroup::signal(std::shared_ptr<Model> m,
                                      const Eigen::VectorXd &p) const {
    Eigen::ArrayXcd result(size());
    size_t start = 0;
    for (auto &sig : sequences) {
        Eigen::ArrayXcd thisResult = sig->signal(m, p);
        result.segment(start, sig->size()) = thisResult;
        start += sig->size();
    }
    return result;
}

Eigen::ArrayXd SequenceGroup::weights(const double f0) const {
    Eigen::ArrayXd weights(size());
    size_t start = 0;
    for (auto &sig : sequences) {
        weights.segment(start, sig->size()) = sig->weights(f0);
        start += sig->size();
    }
    return weights;
}

void SequenceGroup::addSequence(const std::shared_ptr<SequenceBase> &w) {
    sequences.push_back(w);
}

} // End namespace QI
