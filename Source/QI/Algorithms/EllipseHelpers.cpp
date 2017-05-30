/*
 *  EllipseHelpers.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Algorithms/EllipseHelpers.h"

namespace QI {

// Helper Functions
void SemiaxesToHoff(const double A, const double B, const double c,
                    double &a, double &b) {
    b = (-c*A + sqrt(c*c*A*A - (c*c + B*B)*(A*A - B*B)))/(c*c + B*B);
    a = B / (b*B + c*sqrt(1-b*b));
}


void EllipseToMRI(const double a, const double b, const double c, const double th, const double TR, const double flip,
                  float &M0, float &T1, float &T2, float &df0, const bool debug) {
    const double cosf = cos(flip);
    T2 = -TR / log(a);
    T1 = -TR / (log(a-b + a*cosf*(1.-a*b)) - log(a*(1.-a*b) + (a-b)*cosf));
    const double E1 = exp(-TR/T1);
    const double &E2 = a;
    const double M = (c/sqrt(a))*(1-b*b)/(1-a*b);
    M0 = M * (1. - E1*cosf - E2*E2*(E1-cosf)) / ((1-E1)*sin(flip));
    df0 = th / (2*M_PI*TR);
    if (debug) {
        std::cout << "a: " << a << " b: " << b << " c: " << c << " th: " << th << " TR: " << TR << " flip: " << flip << std::endl;
        std::cout << "M0: " << M0 << " T1: " << T1 << " T2: " << T2 << " df0: " << df0 << std::endl;
    }
}

} // End namespace QI