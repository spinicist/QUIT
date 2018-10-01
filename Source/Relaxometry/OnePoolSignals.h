/*
 *  SPGRSignal.h - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_SPGRSIGNAL_H
#define QI_SPGRSIGNAL_H

#include <Eigen/Core>
#include "SPGRSequence.h"
#include "MPRAGESequence.h"
#include "Macro.h"

namespace QI {

// For DESPOT1 B1 is a fixed (double), but for HIFI it is varying (might be a Jet)
template<typename Ta, typename Tb>
inline auto SPGRSignal(const Ta &PD, const Ta &T1, const Tb &B1,
                       const QI::SPGRSequence &s) -> QI_ARRAY(Ta)
{
    const QI_ARRAY(Tb) sa = sin(B1 * s.FA);
    const QI_ARRAY(Tb) ca = cos(B1 * s.FA);
    const Ta E1 = exp(-s.TR / T1);
    return PD * ((1. - E1) * sa) / (1. - E1*ca);
}

template<typename Ta, typename Tb>
inline auto SPGREchoSignal(const Ta &PD, const Ta &T1, const Ta &T2, const Tb &B1,
                           const QI::SPGREchoSequence *s) -> QI_ARRAY(Ta)
{
    const QI_ARRAY(Tb) sa = sin(B1 * s->FA);
    const QI_ARRAY(Tb) ca = cos(B1 * s->FA);
    const Ta E1 = exp(-s->TR / T1);
    return PD * exp(-s->TE / T2) * ((1. - E1) * sa) / (1. - E1*ca);
}

template<typename T> 
inline auto MPRAGESignal(const T &M0, const T &T1, const T &B1,
                         const QI::MPRAGESequence &s) -> QI_ARRAY(T)
{
    const double eta = -1.0; // Inversion efficiency defined as -1 < eta < 0

    const double TIs = s.TI - s.TR*s.k0; // Adjust TI for k0
    const T T1s = 1. / (1./T1 - log(cos(s.FA * B1))/s.TR);
    const T M0s = M0* (1. - exp(-s.TR/T1)) / (1. - exp(-s.TR/T1s));
    const T A_1 = M0s*(1. - exp(-(s.ETL*s.TR)/T1s));

    const T A_2 = M0*(1. - exp(-s.TD/T1));
    const T A_3 = M0*(1. - exp(-TIs/T1));
    const T B_1 = exp(-(s.ETL*s.TR)/T1s);
    const T B_2 = exp(-s.TD/T1);
    const T B_3 = eta*exp(-TIs/T1);

    const T A = A_3 + A_2*B_3 + A_1*B_2*B_3;
    const T B = B_1*B_2*B_3;
    const T M1 = A / (1. - B);

    return QI_ARRAYN(T, 1)(M0s + (M1 - M0s)*exp(-(s.k0*s.TR)/T1s)) * sin(s.FA * B1);
}

} // End namespace QI

#endif // QI_SPGRSIGNAL_H