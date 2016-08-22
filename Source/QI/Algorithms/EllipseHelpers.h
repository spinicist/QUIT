/*
 *  EllipseHelpers.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_ELLIPSE_HELPERS_H
#define QI_ELLIPSE_HELPERS_H

#include <complex>
#include "Eigen/Dense"

#include "QI/Macro.h"
#include "QI/Types.h"
#include "QI/Util.h"

namespace QI {

void SemiaxesToHoff(const double A, const double B, const double c,
                    double &a, double &b);

void EllipseToMRI(const double a, const double b, const double c, const double th, const double TR, const double flip,
                  float &M0, float &T1, float &T2, float &df0);

} // End namespace QI

#endif // QI_ELLIPSEHELPERS_H