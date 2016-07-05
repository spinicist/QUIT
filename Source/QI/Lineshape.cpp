/*
 *  Lineshape.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Lineshape.h"

namespace QI {

Eigen::ArrayXd GaussianLineshape::operator () (const Eigen::ArrayXd &df0, const double T2b) const {
    return sqrt(M_PI_2) * T2b * exp(-pow(2.0*M_PI*df0*T2b,2.0)/2.0);
}

} // End namespace QI
