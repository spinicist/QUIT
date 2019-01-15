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

#ifndef QI_TWOPOOLMODEL_H
#define QI_TWOPOOLMODEL_H

#include "Macro.h"
#include "SPGRSequence.h"
#include "SSFPSequence.h"
#include "SequenceGroup.h"
#include <Eigen/Core>
#include <array>
#include <string>

namespace QI {

struct TwoPoolModel {
    using DataType      = double;
    using ParameterType = double;
    using SequenceType  = QI::SequenceGroup;

    static constexpr int NV = 7;
    static constexpr int ND = 0;
    static constexpr int NF = 2;

    static std::array<const std::string, 7> varying_names;
    static std::array<const std::string, 2> fixed_names;
    static const QI_ARRAYN(double, 2) fixed_defaults;

    const SequenceType &sequence;
    bool                scale_to_mean = false;

    QI_ARRAYN(double, 7) bounds_lo;
    QI_ARRAYN(double, 7) bounds_hi;

    TwoPoolModel(const SequenceType &s, const bool scale);
    bool   valid(const QI_ARRAYN(double, NV) & params) const; // For SRC
    size_t num_outputs() const { return sequence.count(); }
    int    output_size(int i) { return sequence.at(i)->size(); }

    Eigen::MatrixXd SSFP2(const Eigen::ArrayXd &  varying, const QI_ARRAYN(double, NF) & fixed,
                          const QI::SSFPSequence *s) const; // Helper function for SSFP

    Eigen::ArrayXd spgr_signal(const Eigen::ArrayXd &  varying, const QI_ARRAYN(double, NF) & fixed,
                               const QI::SPGRSequence *s) const;

    Eigen::ArrayXd ssfp_signal(const Eigen::ArrayXd &  varying, const QI_ARRAYN(double, NF) & fixed,
                               const QI::SSFPSequence *s) const;

    Eigen::ArrayXd spgr_signal(const Eigen::ArrayXd &varying, const QI_ARRAYN(double, NF) & fixed,
                               const QI::SPGREchoSequence *s) const;

    Eigen::ArrayXd ssfp_signal(const Eigen::ArrayXd &varying, const QI_ARRAYN(double, NF) & fixed,
                               const QI::SSFPEchoSequence *s) const;

    Eigen::ArrayXd signal(const Eigen::ArrayXd &  varying, const QI_ARRAYN(double, NF) & fixed,
                          const QI::SequenceBase *s) const;

    Eigen::ArrayXd signal(const Eigen::ArrayXd &varying, const QI_ARRAYN(double, NF) & fixed) const;

    std::vector<Eigen::ArrayXd> signals(const Eigen::ArrayXd &varying,
                                        const QI_ARRAYN(double, NF) & fixed) const;
};

} // End namespace QI

#endif // QI_TWOPOOLMODEL_H