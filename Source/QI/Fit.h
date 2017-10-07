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
#include <unsupported/Eigen/Splines>

namespace QI {

Eigen::VectorXd LeastSquares(const Eigen::MatrixXd &X, const Eigen::VectorXd &y);
Eigen::VectorXd RobustLeastSquares(const Eigen::MatrixXd &X, const Eigen::VectorXd &y);

/*
 * Eigen spline objects aren't very friendly. Wrap them in a class to do the required
 * scaling and transposes to get them working.
 */
class SplineInterpolator {
public:
    typedef Eigen::Spline<double, 1> TSpline;

    SplineInterpolator() :
        m_min(0.),
        m_width(0.)
    { }

    SplineInterpolator(Eigen::ArrayXd const &x, Eigen::ArrayXd const &y) {
        m_min = x.minCoeff();
        m_width = x.maxCoeff() - m_min;
        const Eigen::ArrayXd sx = (x - m_min) / m_width;
        m_spline = Eigen::SplineFitting<TSpline>::Interpolate(y.transpose(), std::min<int>(x.rows() - 1, 3), sx.transpose());
    }

    double operator()(const double &x) const {
        const double sx = scale(x);
        const double val = m_spline(sx)[0];
        return val;
    }

protected:
    TSpline m_spline;
    double m_min;
    double m_width;

    double scale(const double &x) const {
        return (x - m_min) / m_width;
    }
};


} // End namespace QI

#endif // QI_FIT_H