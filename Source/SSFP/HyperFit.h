/*
 *  HyperFit.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_ELLIPSE_HYPER_H
#define QI_ELLIPSE_HYPER_H

#include "EllipseModel.h"

namespace QI {

struct HyperFit : EllipseFit {
    using EllipseFit::EllipseFit;
    QI::FitReturnType fit(const std::vector<Eigen::ArrayXcd> &inputs,
                          const Eigen::ArrayXd &fixed, QI_ARRAYN(OutputType, EllipseModel::NV) &outputs,
                          ResidualType &residual, std::vector<Eigen::ArrayXcd> &residuals, FlagType &iterations,
                          const int block) const override;
};

} // End namespace QI

#endif // QI_ELLIPSE_HYPER_H