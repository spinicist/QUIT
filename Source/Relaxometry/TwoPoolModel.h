/*
 *  TwoPoolModel.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#pragma once

#include "Model.h"
#include "SPGRSequence.h"
#include "SSFPSequence.h"
#include "SequenceGroup.h"
#include <Eigen/Core>
#include <array>
#include <string>

namespace QI {

struct TwoPoolModel : Model<double, double, 7, 2, 2> {
    std::array<const std::string, 7> const varying_names{
        "PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "tau_m", "f_m"};
    std::array<const std::string, 2> const fixed_names{"f0", "B1"};
    FixedArray const                       fixed_defaults{0.0, 1.0};

    SPGREchoSequence spgr;
    SSFPSequence     ssfp;
    bool             scale_to_mean = false;

    QI_ARRAYN(double, 7) bounds_lo;
    QI_ARRAYN(double, 7) bounds_hi;

    TwoPoolModel(SPGREchoSequence const &s1, SSFPSequence const &s2, const bool scale);
    bool   valid(VaryingArray const &params) const; // For SRC
    size_t num_outputs() const;
    int    output_size(int i) const;

    Eigen::ArrayXd  spgr_signal(const Eigen::ArrayXd &varying,
                                const QI_ARRAYN(double, NF) & fixed) const;
    Eigen::ArrayXd  ssfp_signal(const Eigen::ArrayXd &varying,
                                const QI_ARRAYN(double, NF) & fixed) const;
    Eigen::ArrayXcd ssfp_signalx(const Eigen::ArrayXd &varying,
                                 const QI_ARRAYN(double, NF) & fixed) const;

    std::vector<Eigen::ArrayXd> signals(VaryingArray const &varying, FixedArray const &fixed) const;
    Eigen::ArrayXd signal(const Eigen::ArrayXd &varying, const QI_ARRAYN(double, NF) & fixed) const;
};

} // End namespace QI