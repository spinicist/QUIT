/*
 *  GoldenSection.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "GoldenSection.h"

namespace QI {

double GoldenSectionSearch(std::function<double(double)> f, double a, double b, const double tol) {
    //cout << "goldenSection" << endl;
    const double gr = (std::sqrt(5.0) + 1.0) / 2.0;
    double c = b -  (b - a) / gr;
    double d = a + (b - a) / gr;
    //cout << "c " << c << " d " << d << endl;
    while (std::fabs(c - d) > tol) {
        if (f(c) < f(d)) {
            b = d;
        } else {
            a = c;
        }
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
        //cout << "c " << c << " d " << d << endl;
    }
    return (b + a) / 2.0;
}

} // End namespace QI