/*
 *  SequenceBase.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Sequences/SequenceBase.h"

namespace QI {

/******************************************************************************
 * SequenceBase
 *****************************************************************************/
ostream& operator<<(ostream& os, const SequenceBase& s) {
    s.write(os);
    return os;
}

/******************************************************************************
 * SequenceGroup Class
 *****************************************************************************/
SequenceGroup::SequenceGroup() : SequenceBase()
{}

void SequenceGroup::write(ostream &os) const {
    os << "Combined Sequence Count: " << m_sequences.size() << "\tCombined size: " << size() << endl;
    for (auto& s : m_sequences)
        os << *s;
}

size_t SequenceGroup::count() const {
    return m_sequences.size();
}

shared_ptr<SequenceBase> SequenceGroup::sequence(const size_t i) const {
    return m_sequences.at(i);
}

vector<shared_ptr<SequenceBase>> &SequenceGroup::sequences() {
    return m_sequences;
}

size_t SequenceGroup::size() const {
    size_t sz = 0;
    for (auto& sig : m_sequences)
        sz += sig->size();
    return sz;
}

ArrayXcd SequenceGroup::signal(shared_ptr<Model> m, const VectorXd &p) const {
    ArrayXcd result(size());
    size_t start = 0;
    for (auto &sig : m_sequences) {
        ArrayXcd thisResult = sig->signal(m, p);
        result.segment(start, sig->size()) = thisResult;
        start += sig->size();
    }
    return result;
}

ArrayXd SequenceGroup::weights(const double f0) const {
    ArrayXd weights(size());
    size_t start = 0;
    for (auto &sig : m_sequences) {
        weights.segment(start, sig->size()) = sig->weights(f0);
        start += sig->size();
    }
    return weights;
}

void SequenceGroup::addSequence(const shared_ptr<SequenceBase> &seq) {
    m_sequences.push_back(seq);
}

} // End namespace QI
