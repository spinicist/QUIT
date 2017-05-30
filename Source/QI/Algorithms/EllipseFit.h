/*
 *  Ellipse.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_ELLIPSEFIT_H
#define QI_ELLIPSEFIT_H

#include "QI/Algorithms/Ellipse.h"

namespace QI {

class FitEllipse : public ESAlgo {
public:
    FitEllipse(std::shared_ptr<QI::SSFPEcho> &seq, bool debug, bool phase, bool block) :
        ESAlgo(seq, debug, phase), m_reorderBlock(block)
    {}

protected:
    bool m_reorderBlock;

    virtual std::array<float, NumOutputs> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &indata,
                                           const double TR, const double flip) const;
};

} // End namespace QI

#endif // QI_ELLIPSEFIT_H