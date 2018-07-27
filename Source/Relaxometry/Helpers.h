/*
 *  Helpers.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_RELAX_HELPERS_H
#define QI_RELAX_HELPERS_H

#include <Eigen/Core>

namespace Eigen {
    using Matrix6d = Matrix<double, 6, 6>;
    using Vector6d = Matrix<double, 6, 1>;
    using Matrix9d = Matrix<double, 9, 9>;
    using Vector9d = Matrix<double, 9, 1>;
    using MagVector = Matrix<double, 3, Eigen::Dynamic>;
}

namespace QI {

void CalcExchange(const double tau_a, const double f_a, double &f_b, double &k_ab, double &k_ba);

} // End namespace QI

#endif // QI_RELAX_HELPERS_H
