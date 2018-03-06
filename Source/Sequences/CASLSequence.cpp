/*
 *  CASLSequence.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "CASLSequence.h"
#include "Macro.h"
#include "EigenCereal.h"

namespace QI {

void CASLSequence::load(cereal::JSONInputArchive &ar) {
    ar(CEREAL_NVP(TR),
       CEREAL_NVP(label_time),
       CEREAL_NVP(post_label_delay));
}

void CASLSequence::save(cereal::JSONOutputArchive &ar) const {
    ar(CEREAL_NVP(TR),
       CEREAL_NVP(label_time),
       CEREAL_NVP(post_label_delay));
}

Eigen::ArrayXcd CASLSequence::signal(const std::shared_ptr<QI::Model> m, const Eigen::VectorXd &p) const {
    QI_FAIL("Not Implemented");
}

} // End namespace QI
