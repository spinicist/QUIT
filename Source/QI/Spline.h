/*
 *  Spline.h
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, you can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_SPLINE_H
#define QI_SPLINE_H

#include <iostream>

#include "Eigen/Dense"
#include <unsupported/Eigen/Splines>

namespace QI {

/*
 * Eigen spline objects aren't very friendly. Wrap them in a class to do the required
 * scaling and transposes to get them working.
 */
class SplineInterpolator {
public:
    typedef Eigen::Spline<double, 1> TSpline;

    SplineInterpolator();
    SplineInterpolator(Eigen::ArrayXd const &x, Eigen::ArrayXd const &y);
    double operator()(const double &x) const;

protected:
    TSpline m_spline;
    double m_min;
    double m_width;
    double scale(const double &x) const;
};

} // End namespace QI

#endif // QI_SPLINE_H