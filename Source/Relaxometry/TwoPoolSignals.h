/*
 *  TwoPoolSignals.h - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2021 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#pragma once

#include "Macro.h"
#include "SPGRSequence.h"
#include "SSFPSequence.h"
#include <Eigen/Core>

namespace QI {

Eigen::ArrayXcd SPGR2(double const            PD,
                      double const            T1_a,
                      double const            T2_a,
                      double const            T1_b,
                      double const            T2_b,
                      double const            tau_a,
                      double const            f_a,
                      double const            f0,
                      double const            B1,
                      SPGREchoSequence const &spgr);

Eigen::ArrayXcd SSFP2(double const        PD,
                      double const        T1_a,
                      double const        T2_a,
                      double const        T1_b,
                      double const        T2_b,
                      double const        tau_a,
                      double const        f_a,
                      double const        f0,
                      double const        B1,
                      SSFPSequence const &ssfp);

} // namespace QI
