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

#include "SequenceBase.h"

namespace QI {

size_t SequenceBase::count() const{
    return 1;
}

Eigen::ArrayXd SequenceBase::weights(const double /* Unused */) const {
    return Eigen::ArrayXd::Ones(size()); // Default weights are constant
}

Eigen::ArrayXd SequenceBase::signal_magnitude(const std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
    Eigen::ArrayXcd c_signal = this->signal(m, p);
    return c_signal.abs();
}

} // End namespace QI
