/*
 *  ThreePoolModel.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_THREEPOOLMODEL_H
#define QI_THREEPOOLMODEL_H

#include "Model.h"
#include "OnePoolSignals.h"
#include "SPGRSequence.h"
#include "SSFPSequence.h"
#include <Eigen/Core>
#include <array>
#include <string>

namespace QI {

struct ThreePoolModel : Model<double, double, 10, 2, 2> {
    SPGREchoSequence spgr;
    SSFPSequence     ssfp;
    bool             scale_to_mean = false;

    std::array<const std::string, ThreePoolModel::NV> const varying_names{
        "PD", "T1_m", "T2_m", "T1_ie", "T2_ie", "T1_csf", "T2_csf", "tau_m", "f_m", "f_csf"};
    std::array<const std::string, ThreePoolModel::NF> const fixed_names{"f0", "B1"};
    FixedArray const                                        fixed_defaults{0.0, 1.0};

    QI_ARRAYN(double, NV) bounds_lo;
    QI_ARRAYN(double, NV) bounds_hi;

    ThreePoolModel(SPGREchoSequence const &s1, SSFPSequence const &s2, const bool scale);
    bool   valid(const QI_ARRAYN(double, NV) & params) const; // For SRC
    size_t num_outputs() const;
    int    output_size(int i) const;

    Eigen::ArrayXd spgr_signal(const Eigen::ArrayXd &varying,
                               const QI_ARRAYN(double, NF) & fixed) const;

    Eigen::ArrayXd ssfp_signal(const Eigen::ArrayXd &varying,
                               const QI_ARRAYN(double, NF) & fixed) const;

    std::vector<Eigen::ArrayXd> signals(const Eigen::ArrayXd &varying,
                                        const QI_ARRAYN(double, NF) & fixed) const;
    Eigen::ArrayXd signal(const Eigen::ArrayXd &varying, const QI_ARRAYN(double, NF) & fixed) const;
};

} // End namespace QI

#endif // QI_THREEPOOLMODEL_H