/*
 *  Helpers.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <limits>
#include <cmath>

namespace QI {

// Helper Functions
// Calculate the exchange rates from the residence time and fractions
void CalcExchange(const double tau_a, const double f_a, double &f_b, double &k_ab, double &k_ba) {
    const double feps = std::numeric_limits<float>::epsilon(); // Because we read from float files
    f_b = 1.0 - f_a;
    k_ab = 1./tau_a; k_ba = k_ab*f_a/f_b;
    if ((std::fabs(f_a - 1.) <= feps) || (std::fabs(f_b - 1.) <= feps)) {
        // Only have 1 component, so no exchange
        k_ab = 0.;
        k_ba = 0.;
    }
}

} // End namespace QI
