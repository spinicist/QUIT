/*
 *  TwoPoolModel.hxx
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  Based on code by Sean Deoni
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */
// #define QI_DEBUG_BUILD
#include "Macro.h"

#include "Helpers.h"
#include "Log.h"
#include "TwoPoolModel.h"
#include "TwoPoolSignals.h"

// There's a bug in matrix log that leads to duplicate symbol definitions if
// this is included in the header file
#include <unsupported/Eigen/MatrixFunctions>

using namespace std::literals;

namespace QI {

TwoPoolModel::TwoPoolModel(const SPGREchoSequence &s1, const SSFPSequence &s2, const bool scale) :
    spgr{s1}, ssfp{s2}, scale_to_mean{scale} {
    bounds_lo << 1.0, 0.300, 0.010, 0.9, 0.040, 0.025, 0.001;
    bounds_hi << 1.0, 0.800, 0.030, 3.0, 1.500, 0.600, 0.35;
}

size_t TwoPoolModel::num_outputs() const {
    return 2;
}

int TwoPoolModel::output_size(int i) const {
    if (i == 0) {
        return spgr.size();
    } else if (i == 1) {
        return ssfp.size();
    } else {
        QI::Fail("Invalid output size: {}", i);
    }
}

bool TwoPoolModel::valid(VaryingArray const &params) const {
    // Negative T1/T2 makes no sense
    if ((params[1] <= 0.) || (params[2] <= 0.))
        return false;
    else {
        if ((params[1] < params[3]) && (params[2] < params[4]) && (params[6] <= 1.0))
            return true;
        else
            return false;
    }
}

std::vector<Eigen::ArrayXd> TwoPoolModel::signals(VaryingArray const &v,
                                                  FixedArray const &  f) const {
    return {spgr_signal(v, f), ssfp_signal(v, f)};
}

Eigen::ArrayXd TwoPoolModel::signal(const Eigen::ArrayXd &v,
                                    const QI_ARRAYN(double, NF) & f) const {
    auto           sigs = signals(v, f);
    Eigen::ArrayXd sig(spgr.size() + ssfp.size());
    sig.head(spgr.size()) = sigs[0];
    sig.tail(ssfp.size()) = sigs[1];
    return sig;
}

Eigen::ArrayXd TwoPoolModel::spgr_signal(const Eigen::ArrayXd &v,
                                         const QI_ARRAYN(double, NF) & f) const {
    Eigen::ArrayXcd const M = SPGR2(v[0], v[1], v[2], v[3], v[4], v[5], v[6], f[0], f[1], spgr);
    Eigen::ArrayXd        signal = M.abs();
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}

Eigen::ArrayXd TwoPoolModel::ssfp_signal(const Eigen::ArrayXd &v,
                                         const QI_ARRAYN(double, NF) & f) const {
    Eigen::ArrayXcd const M = SSFP2(v[0], v[1], v[2], v[3], v[4], v[5], v[6], f[0], f[1], ssfp);
    Eigen::ArrayXd        signal = M.abs();
    QI_DBMAT(M);
    QI_DBVEC(signal);
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}

} // End namespace QI
