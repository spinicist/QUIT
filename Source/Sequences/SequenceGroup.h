/*
 *  SequenceGroup.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_GROUP_H
#define SEQUENCES_GROUP_H

#include "Macro.h"
#include "SequenceBase.h"
#include <vector>

namespace QI {

struct SequenceGroup : SequenceBase {
    std::vector<const SequenceBase *> sequences;
    const SequenceBase *              at(const size_t i) const;
    const SequenceBase *              operator[](const size_t i) const;

    std::string &name() const override {
        static std::string name = "SequenceGroup";
        return name;
    }
    size_t         count() const override;
    Eigen::Index   size() const override;
    Eigen::ArrayXd weights(const double f0 = 0.0) const override;

    void addSequence(const SequenceBase *s);
};

} // End namespace QI

#endif // SEQUENCES_BASE_H
