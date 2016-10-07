/*
 *  CEST.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Algorithms/CEST.h"

namespace QI {

void CESTAlgo::apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                     std::vector<TOutput> &outputs, TConst &residual,
                     TInput &resids, TIters &its) const
{
    size_t full = m_half*2+1;
    const Eigen::Map<const Eigen::ArrayXf> z_spec(inputs[0].GetDataPointer(), full);
    for (int f = 0; f < m_half; f++) {
        outputs[0][f] = (z_spec[m_half + f] - z_spec[m_half - 1 - f]) / z_spec[0];
    }
}

} // End namespace QI