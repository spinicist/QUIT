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

#ifndef QI_ELLIPSEFITS_H
#define QI_ELLIPSEFITS_H

#include <Eigen/Dense>

namespace QI {

typedef Eigen::Array<double, 5, 1> Array5d;

enum class EllipseMethods { Hyper, Direct };

Array5d HyperEllipse(const Eigen::ArrayXcf &input, const double TR, const Eigen::ArrayXd &phi);
Array5d DirectEllipse(const Eigen::ArrayXcf &input, const double TR, const Eigen::ArrayXd &phi, const bool debug);

} // End namespace QI

#endif // QI_ELLIPSEFIT_H