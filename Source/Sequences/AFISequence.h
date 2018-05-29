/*
 *  AFI.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_AFI_H
#define SEQUENCES_AFI_H

#include "SequenceBase.h"
#include <cereal/cereal.hpp>

namespace QI {

struct AFISequence : SequenceBase {
    double FA, TR1, TR2;
    Eigen::Index size() const override { return 2; }
    Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;

    template<typename Archive>
    void serialize(Archive &archive) {
        archive(CEREAL_NVP(FA), CEREAL_NVP(TR1), CEREAL_NVP(TR2));
    }
};

} // End namespace QI

#endif // SEQUENCES_AFI_H