/*
 *  OnePoolSignals.h - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2021 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#pragma once

#include "MPRAGESequence.h"
#include "Macro.h"
#include "SPGRSequence.h"
#include "SSFPSequence.h"
#include <Eigen/Core>

namespace QI {

Eigen::ArrayXcd SPGR1(double const &          PD,
                      double const &          T1,
                      double const &          T2,
                      double const &          f0,
                      double const &          B1,
                      SPGREchoSequence const &s);
Eigen::ArrayXcd SSFP1(double const &          PD,
                      double const &          T1,
                      double const &          T2,
                      double const &          f0,
                      double const &          B1,
                      const QI::SSFPSequence &s);

// For DESPOT1 B1 is a fixed (double), but for HIFI it is varying (might be a Jet)
template <typename Ta, typename Tb>
inline auto SPGRSignal(Ta const &PD, Ta const &T1, Tb const &B1, SPGRSequence const &s)
    -> QI_ARRAY(Ta) {
    const QI_ARRAY(Tb) sa = sin(B1 * s.FA);
    const QI_ARRAY(Tb) ca = cos(B1 * s.FA);
    Ta const E1           = exp(-s.TR / T1);
    return PD * ((1. - E1) * sa) / (1. - E1 * ca);
}

template <typename T>
inline auto MPRAGESignal(T const &M0, T const &T1, T const &B1, MPRAGESequence const &s)
    -> QI_ARRAY(T) {
    double const eta = -1.0; // Inversion efficiency defined as -1 < eta < 0

    double const TIs = s.TI - s.TR * s.k0; // Adjust TI for k0
    T const      T1s = 1. / (1. / T1 - log(cos(s.FA * B1)) / s.TR);
    T const      M0s = M0 * (1. - exp(-s.TR / T1)) / (1. - exp(-s.TR / T1s));
    T const      A_1 = M0s * (1. - exp(-(s.ETL * s.TR) / T1s));

    T const A_2 = M0 * (1. - exp(-s.TD / T1));
    T const A_3 = M0 * (1. - exp(-TIs / T1));
    T const B_1 = exp(-(s.ETL * s.TR) / T1s);
    T const B_2 = exp(-s.TD / T1);
    T const B_3 = eta * exp(-TIs / T1);

    T const A  = A_3 + A_2 * B_3 + A_1 * B_2 * B_3;
    T const B  = B_1 * B_2 * B_3;
    T const M1 = A / (1. - B);

    return QI_ARRAYN(T, 1)(M0s + (M1 - M0s) * exp(-(s.k0 * s.TR) / T1s)) * sin(s.FA * B1);
}

} // End namespace QI
