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
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(X.rows(), X.cols());
    QR.compute(X);
    Eigen::VectorXd b = QR.solve(y);
    Eigen::ArrayXd corr_fac;
    {   // Calculate leverage and correction factor from the 'thinQ' matrix in here to de-alloc memory when done
        Eigen::MatrixXd Q(X.rows(), X.cols()); Q.setIdentity();
        Q = QR.householderQ() * Q;
        // The below expression with Q is equal to (Q * Q.transpose()).diagonal() = h
        corr_fac = (1.0 - Q.array().square().rowwise().sum()).sqrt();
    }

    // Allocate some workspace for loop
    Eigen::ArrayXd r(y.rows()), sr(y.rows()), w(y.rows());
    Eigen::VectorXd b_prev(b.rows());
    Eigen::MatrixXd wX(X.rows(), X.cols());
    
    const double tune = 1.345; // For Huber only
    int iter = 0;
    bool converged = false;
    while (!converged && (++iter < 20)) {
        r = (y - X*b); // Calculate residuals
        const double sig = mad_sigma(r, sr, X.cols()); // Get Median Absolute Deviation
        r /= (tune * std::max(sig, sig_lower) * corr_fac); // Adjust residuals
        w = 1 / r.abs().max(1); // Huber weights
        b_prev = b; // Save co-efficients
        wX = w.matrix().asDiagonal() * X;
        QR.compute(wX);
        w *= y.array(); // Re-use w array for weighted-y
        b = QR.solve(w.matrix()); // Weighted solve
        if (((b - b_prev).array().abs() < sqrt(std::numeric_limits<double>::epsilon())).all()) {
            converged = true;
        }
    }

    return b;
}

} // End namespace QI