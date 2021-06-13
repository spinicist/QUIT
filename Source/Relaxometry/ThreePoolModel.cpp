/*
 *  ThreePoolModel.hxx
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

#include "Helpers.h"
#include "Log.h"
#include "OnePoolSignals.h"
#include "ThreePoolModel.h"
#include "TwoPoolSignals.h"

#include <Eigen/Dense>
// There's a bug in matrix log that leads to duplicate symbol definitions if
// this is included in the header file
#include <unsupported/Eigen/MatrixFunctions>

using namespace std::literals;

namespace QI {

/* The inner two-pool model must not be scaled, as scaling will be done once signals are added */
ThreePoolModel::ThreePoolModel(SPGREchoSequence const &s1,
                               SSFPSequence const &    s2,
                               const bool              scale) :
    spgr{s1},
    ssfp{s2}, scale_to_mean{scale} {
    bounds_lo << 1.0, 0.300, 0.010, 0.9, 0.040, 3.5, 1.0, 0.025, 0.001, 0.001;
    bounds_hi << 1.0, 0.800, 0.030, 1.5, 0.150, 5.0, 3.5, 0.600, 0.350, 0.999;
}

size_t ThreePoolModel::num_outputs() const {
    return 2;
}

int ThreePoolModel::output_size(int i) const {
    if (i == 0) {
        return spgr.size();
    } else if (i == 1) {
        return ssfp.size();
    } else {
        QI::Fail("Invalid output size: {}", i);
    }
}

bool ThreePoolModel::valid(const QI_ARRAYN(double, NV) & v) const {
    // Negative T1/T2 makes no sense
    if ((v[1] <= 0.) || (v[2] <= 0.))
        return false;
    else {
        if ((v[1] < v[3]) && (v[2] < v[4]) && (v[3] < v[5]) && (v[4] < v[6]) &&
            ((v[8] + v[9]) <= 1.0))
            return true;
        else
            return false;
    }
}

std::vector<Eigen::ArrayXd> ThreePoolModel::signals(const Eigen::ArrayXd &v,
                                                    const QI_ARRAYN(double, NF) & f) const {
    return {spgr_signal(v, f), ssfp_signal(v, f)};
}

Eigen::ArrayXd ThreePoolModel::signal(const Eigen::ArrayXd &v,
                                      const QI_ARRAYN(double, NF) & f) const {
    auto           sigs   = signals(v, f);
    Eigen::ArrayXd sig    = Eigen::ArrayXd::Zero(spgr.size() + ssfp.size());
    sig.head(spgr.size()) = sigs[0];
    sig.tail(ssfp.size()) = sigs[1];
    QI_DBVEC(sig);
    return sig;
}

Eigen::ArrayXd ThreePoolModel::spgr_signal(const Eigen::ArrayXd &v,
                                           const QI_ARRAYN(double, NF) & fixed) const {
    double     f_ab = 1. - v[9];
    auto const m_ab =
        SPGR2(v[0] * f_ab, v[1], v[2], v[3], v[4], v[7], v[8] / f_ab, fixed[0], fixed[1], spgr);
    auto const     m_c    = SPGR1(v[0] * v[9], v[5], v[6], fixed[0], fixed[1], spgr);
    Eigen::ArrayXd signal = (m_ab + m_c).abs();
    QI_DBMSG("spgr\n");
    QI_DBVEC(v);
    QI_DBVEC(m_ab);
    QI_DBVEC(m_c);
    QI_DBVEC(signal);
    // std::cout << "SPGR\n" << m_ab.transpose() << "\n" << m_c.transpose() << "\n" <<
    // signal.transpose() << std::endl;
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}

Eigen::ArrayXd ThreePoolModel::ssfp_signal(const Eigen::ArrayXd &v,
                                           const QI_ARRAYN(double, NF) & fixed) const {
    double     f_ab = 1. - v[9];
    auto const m_ab =
        SSFP2(v[0] * f_ab, v[1], v[2], v[3], v[4], v[7], v[8] / f_ab, fixed[0], fixed[1], ssfp);
    auto const     m_c    = SSFP1(v[0] * v[9], v[5], v[6], fixed[0], fixed[1], ssfp);
    Eigen::ArrayXd signal = (m_ab + m_c).abs();
    QI_DBMSG("ssfp\n");
    QI_DBVEC(two_pool_varying);
    QI_DBVEC(v);
    QI_DBVEC(m_ab);
    QI_DBVEC(m_c);
    QI_DBVEC(signal);
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}

} // End namespace QI
