/*
 *  HyperAlgo.h
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_ELLIPSE_HYPER_H
#define QI_ELLIPSE_HYPER_H

#include <Eigen/Dense>
#include "EllipseAlgo.h"

namespace QI {

class HyperAlgo : public EllipseAlgo {
protected:
    Eigen::ArrayXd apply_internal(const Eigen::ArrayXcf &input, const double flip, const double TR, const Eigen::ArrayXd &phi, const bool debug, float &residual) const override;
public:
    HyperAlgo(const QI::SSFPEllipseSequence &seq, bool debug) : EllipseAlgo(seq, debug) {};
};

} // End namespace QI

#endif // QI_ELLIPSE_HYPER_H