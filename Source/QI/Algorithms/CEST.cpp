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

#include <unsupported/Eigen/Splines>
#include "QI/Algorithms/CEST.h"

namespace QI {

void CESTAlgo::apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                     std::vector<TOutput> &outputs, TConst &residual,
                     TInput &resids, TIters &its) const
{
    size_t full = m_half*2+1;
    const Eigen::Map<const Eigen::ArrayXf> z_spec(inputs[0].GetDataPointer(), full);
    size_t index; z_spec.minCoeff(&index);
    float f0 = m_ifrqs.coeffRef(index);
    outputs.at(2)[0] = f0;

    typedef Eigen::Spline<float, 1> TSpline;
    typedef Eigen::SplineFitting<TSpline> TFit;
    const float maxfrq = m_ifrqs.maxCoeff();
    const float minfrq = m_ifrqs.minCoeff();
    const float w = maxfrq - minfrq;
    const Eigen::ArrayXf scaledfrqs = (m_ifrqs - minfrq) / w;
    Eigen::DenseIndex degree = std::min<int>(m_ifrqs.rows() - 1, 3);
    TSpline spline = TFit::Interpolate(z_spec.transpose(), degree, scaledfrqs);
    const float ref = spline(0)[0];
    for (int f = 0; f < m_ifrqs.rows(); f++) {
        const float frq = (m_ifrqs.coeffRef(f) + f0 - minfrq)/w;
        const TSpline::PointType val = spline(frq);
        outputs.at(0)[f] = val[0];
    }
    for (int f = 0; f < m_ofrqs.rows(); f++) {
        const float pfrq = (f0 + m_ofrqs.coeffRef(f) - minfrq)/w;
        const float nfrq = (f0 - m_ofrqs.coeffRef(f) - minfrq)/w;
        const float pos = spline(pfrq)[0];
        const float neg = spline(nfrq)[0];
        outputs.at(1)[f] = ((pos - neg)/ref);
    }
}

} // End namespace QI