/*
 *  GridInterp.h
 *
 *  Copyright (c) 2025 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, you can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_GRIDINTERP_H
#define QI_GRIDINTERP_H

#include "Eigen/Core"

namespace QI {

struct RegularGrid {
    RegularGrid(Eigen::ArrayXd const &x, Eigen::ArrayXd const &y, Eigen::ArrayXXd const &z);
    double operator()(const double &x, const double &y) const;

  protected:
    struct IndexPair { Eigen::Index lo, hi; };
    static auto FindIndices(Eigen::ArrayXd const &a, double const x) -> IndexPair;
    Eigen::ArrayXd  m_x, m_y;
    Eigen::ArrayXXd m_z;
    Eigen::Index const nX, nY;
};

} // End namespace QI

#endif // QI_SPLINE_H