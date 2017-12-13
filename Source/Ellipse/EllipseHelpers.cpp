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

#include "EllipseHelpers.h"

namespace QI {

// Helper Functions

Eigen::ArrayXd Unwrap(const Eigen::ArrayXd &x) {
    // if (debug) {
    //     std::cout << "*** START ***\nY wrapped:   " << Y.transpose() << std::endl;
    // }
    const int sz = x.rows();
    Eigen::ArrayXd diff = x.tail(sz - 1) - x.head(sz - 1);
    const Eigen::ArrayXd diffmod = ((diff+M_PI) - (2*M_PI)*floor((diff+M_PI)/(2*M_PI))) - M_PI;
    Eigen::ArrayXd ph_correct = diffmod - diff;
    for (int i = 0; i < (sz - 1); i++) {
        if (std::abs(diff[i]) < M_PI) ph_correct[i] = 0;
    }
    for (int i = 1; i < (sz - 1); i++) {
        ph_correct[i] += ph_correct[i - 1];
    }
    Eigen::ArrayXd y(x.rows());
    y[0] = x[0];
    y.tail(sz - 1) += ph_correct;
    return y;
    // if (debug) {
    //     std::cout << "Diff         " << diff.transpose() << std::endl;
    //     std::cout << "Diffmod      " << diffmod.transpose() << std::endl;
    //     std::cout << "Ph correct   " << ph_correct.transpose() << std::endl;
    //     std::cout << "Y unwrapped: " << Y.transpose() << std::endl;
    // }
}


void SemiaxesToHoff(const double A, const double B, const double c,
                    double &G, double &a, double &b) {
    b = (-c*A + sqrt(c*c*A*A - (c*c + B*B)*(A*A - B*B)))/(c*c + B*B);
    a = B / (b*B + c*sqrt(1-b*b));
    G = c*(1 - b*b)/(1 - a*b);
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