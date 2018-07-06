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

#include "SequenceBase.h"
#include "Macro.h"

namespace QI {

struct SequenceGroup : SequenceBase {
    std::vector<std::shared_ptr<SequenceBase>> sequences;
    std::shared_ptr<SequenceBase> &operator[](const size_t i);

    std::string &name() const override { static std::string name = "SequenceGroup"; return name; }
    size_t count() const override;
    Eigen::Index size() const override;
    Eigen::ArrayXcd signal(std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &par) const override;
    Eigen::ArrayXd weights(const double f0 = 0.0) const override;

    void addSequence(const std::shared_ptr<SequenceBase> &s);
    SequenceGroup(rapidjson::Value &);
    rapidjson::Value toJSON(rapidjson::Document::AllocatorType &) const override;
};

} // End namespace QI

#endif // SEQUENCES_BASE_H
