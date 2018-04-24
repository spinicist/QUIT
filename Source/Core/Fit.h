/*
 *  Fit.h
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, you can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_FIT_H
#define QI_FIT_H

#include <iostream>

#include "Eigen/Dense"

namespace QI {

Eigen::VectorXd LeastSquares(const Eigen::MatrixXd &X, const Eigen::VectorXd &y, double *resid = nullptr);
Eigen::VectorXd RobustLeastSquares(const Eigen::MatrixXd &X, const Eigen::VectorXd &y, double *resid = nullptr);

} // End namespace QI

#endif // QI_FIT_H