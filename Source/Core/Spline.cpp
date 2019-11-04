/*
 *  Spline.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, you can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Log.h"
#include "Spline.h"

#include <set>

namespace QI {

SplineInterpolator::SplineInterpolator() : m_min(0.), m_width(0.) {}

// Sort input in ascending order

SplineInterpolator::SplineInterpolator(Eigen::ArrayXd const &     x,
                                       Eigen::ArrayXd const &     y,
                                       const int                  order,
                                       std::vector<size_t> const &indices) {
    if (x.size() != y.size()) {
        QI::Fail("Input vectors to spline must be same size");
    }
    if (x.size() == 0) {
        QI::Fail("Cannot create a spline with no control points");
    }

    Eigen::ArrayXd sx, sy;
    if (indices.size() > 0) {
        sx = Eigen::ArrayXd(indices.size());
        sy = Eigen::ArrayXd(indices.size());
        std::transform(indices.begin(), indices.end(), &sx[0], [&](std::size_t i) { return x[i]; });
        std::transform(indices.begin(), indices.end(), &sy[0], [&](std::size_t i) { return y[i]; });
    } else {
        sx = x;
        sy = y;
    }

    // fmt::print("Sorted freqs: {}\nSorted value: {}\n", sx.transpose(), sy.transpose());

    m_min                        = sx[0];
    m_width                      = sx[sx.size() - 1] - m_min;
    const Eigen::ArrayXd scaledx = (sx - m_min) / m_width;
    m_spline                     = Eigen::SplineFitting<TSpline>::Interpolate(
        sy.transpose(), std::min<int>(x.rows() - 1, order), scaledx.transpose());
}

double SplineInterpolator::scale(const double &x) const {
    return (x - m_min) / m_width;
}

Eigen::ArrayXd SplineInterpolator::scale(const Eigen::ArrayXd &x) const {
    return (x - m_min) / m_width;
}

double SplineInterpolator::operator()(const double &x) const {
    const double sx  = scale(x);
    const double val = m_spline(sx)[0];
    return val;
}

Eigen::ArrayXd SplineInterpolator::operator()(const Eigen::ArrayXd &x) const {
    const auto     sx = scale(x);
    Eigen::ArrayXd output(sx.rows());
    for (Eigen::Index i = 0; i < sx.rows(); i++) {
        output[i] = m_spline(sx[i])[0];
    }
    return output;
}

void SplineInterpolator::print(std::ostream &ostr) const {
    ostr << "SPLINE Min: " << m_min << " Width: " << m_width << std::endl;
}

std::ostream &operator<<(std::ostream &ostr, const SplineInterpolator &sp) {
    sp.print(ostr);
    return ostr;
}

} // End namespace QI