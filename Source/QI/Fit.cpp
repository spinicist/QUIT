/*
 *  Fit.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, you can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Fit.h"

#include <iostream>
#include <iomanip>
#include <algorithm>

namespace QI {

double standard_dev(const Eigen::ArrayXd &x) {
    return sqrt((x - x.mean()).square().sum() / (x.rows() - 1));
}

double mad_sigma(Eigen::ArrayXd &r, Eigen::ArrayXd &sr, const int p) {
    sr = r.abs();
    std::sort(sr.data(), sr.data() + sr.size());
    const size_t index = (sr.size() + p) / 2;
    return sr[index] / 0.6745;
}

Eigen::VectorXd LeastSquares(const Eigen::MatrixXd &X, const Eigen::VectorXd &y) {
    // Solve Xb = Y via least squares via QR decomposition
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(X);
    Eigen::VectorXd b = QR.solve(y);
    return b;
}

Eigen::VectorXd RobustLeastSquares(const Eigen::MatrixXd &X, const Eigen::VectorXd &y) {
    // With thanks to gsl_multifit_robust & Matlab
    const double sig_y = standard_dev(y);
    const double sig_lower = (sig_y == 0) ? 1.0 : 1e-6 * sig_y;

    // Solve Xb = Y via least squares via QR decomposition, keep the decomposition around
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(X);
    Eigen::VectorXd b = QR.solve(y);
    // Calculate leverage and correction factor from the 'thinQ' matrix
    Eigen::MatrixXd Q(X.rows(), X.cols()); Q.setIdentity();
    Q = QR.householderQ() * Q;
    Eigen::ArrayXd h = Q.array().square().rowwise().sum(); //(Q * Q.transpose()).diagonal();
    Eigen::ArrayXd corr_fac = (1.0 - h).sqrt();

    // For Huber at the moment
    const double tune = 1.345;
    int iter = 0;
    bool converged = false;

    // Allocate some workspace for loop
    Eigen::ArrayXd resid(y.rows()), sorted_resid(y.rows()), r(y.rows()), w(y.rows()), wy(y.rows());
    Eigen::VectorXd b_prev(b.rows());
    Eigen::MatrixXd wX(X.rows(), X.cols());
    while (!converged && (++iter < 10)) {
        resid = (y - X*b); // Calculate residuals
        const double sig = mad_sigma(resid, sorted_resid, X.cols()); // Get Median Absolute Deviation
        r = resid / (tune * std::max(sig, sig_lower) * corr_fac); // Adjust residuals
        std::cerr << "Iteration: " << iter << " total resid: " << resid.square().sum() << " adjusted resid: " << r.square().sum() << std::endl;
        w = 1 / r.abs().max(1); // Huber weights
        b_prev = b; // Save weights
        wX = w.matrix().asDiagonal() * X;
        wy = w * y.array();
        b = LeastSquares(wX, wy); // Weighted solve
        if (((b - b_prev).array().abs() < sqrt(std::numeric_limits<double>::epsilon())).all()) {
            converged = true;
        }
    }

    return b;
}

} // End namespace QI