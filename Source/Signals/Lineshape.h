/*
 *  Lineshape.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef LINESHAPE_H
#define LINESHAPE_H

#include <functional>
#include <Eigen/Dense>

namespace QI {

typedef std::function<Eigen::ArrayXd(const Eigen::ArrayXd &, const double)> TLineshape;

class GaussianLineshape {
public:
    Eigen::ArrayXd operator()(const Eigen::ArrayXd &df0, const double T2b) const;
};

} // End namespace QI

#endif // LINESHAPE_H
