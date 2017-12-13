/*
 *  GoldenSection.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_GOLDENSECTION_H
#define QI_GOLDENSECTION_H

#include <cmath>
#include <functional>

namespace QI {

double GoldenSectionSearch(std::function<double(double)> f, double a, double b, const double tol);

} // End namespace QI

#endif // QI_GOLDENSECTION_H