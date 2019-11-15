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

#include "Helpers.h"
#include "Log.h"
#include "TwoPoolModel.h"

// There's a bug in matrix log that leads to duplicate symbol definitions if
// this is included in the header file
#include <unsupported/Eigen/MatrixFunctions>

using namespace std::literals;

namespace {
Eigen::MatrixXd SSFP2(const Eigen::ArrayXd &varying,
                      const QI_ARRAYN(double, 2) & fixed,
                      QI::SSFPSequence const &ssfp) {
    const double &PD    = varying[0];
    const double &T1_a  = varying[1];
    const double &T2_a  = varying[2];
    const double &T1_b  = varying[3];
    const double &T2_b  = varying[4];
    const double &tau_a = varying[5];
    const double &f_a   = varying[6];
    const double &f0    = fixed[0];
    const double &B1    = fixed[1];
    const double &TR    = ssfp.TR;
    const double  E1_a  = exp(-TR / T1_a);
    const double  E1_b  = exp(-TR / T1_b);
    const double  E2_a  = exp(-TR / T2_a);
    const double  E2_b  = exp(-TR / T2_b);
    double        f_b, k_ab, k_ba;
    QI::CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    const double         E_ab  = exp(-TR * k_ab / f_b);
    const double         K1    = E_ab * f_b + f_a;
    const double         K2    = E_ab * f_a + f_b;
    const double         K3    = f_a * (1 - E_ab);
    const double         K4    = f_b * (1 - E_ab);
    const Eigen::ArrayXd alpha = B1 * ssfp.FA;
    const Eigen::ArrayXd theta = ssfp.PhaseInc + 2. * M_PI * f0 * TR;

    Eigen::MatrixXd M(4, ssfp.size());
    Eigen::Matrix6d LHS;
    Eigen::Vector6d RHS;
    RHS << 0, 0, 0, 0, -E1_b * K3 * f_b + f_a * (-E1_a * K1 + 1),
        -E1_a * K4 * f_a + f_b * (-E1_b * K2 + 1);

    for (int i = 0; i < ssfp.size(); i++) {
        const double ca  = cos(alpha[i]);
        const double sa  = sin(alpha[i]);
        const double cta = cos(theta[i]);
        const double ctb = cos(theta[i]);
        const double sta = sin(theta[i]);
        const double stb = sin(theta[i]);

        LHS << -E2_a * K1 * cta + ca, -E2_b * K3 * cta, E2_a * K1 * sta, E2_b * K3 * sta, sa, 0,
            -E2_a * K4 * ctb, -E2_b * K2 * ctb + ca, E2_a * K4 * stb, E2_b * K2 * stb, 0, sa,
            -E2_a * K1 * sta, -E2_b * K3 * sta, -E2_a * K1 * cta + 1, -E2_b * K3 * cta, 0, 0,
            -E2_a * K4 * stb, -E2_b * K2 * stb, -E2_a * K4 * ctb, -E2_b * K2 * ctb + 1, 0, 0, -sa,
            0, 0, 0, -E1_a * K1 + ca, -E1_b * K3, 0, -sa, 0, 0, -E1_a * K4, -E1_b * K2 + ca;
        M.col(i).noalias() = PD * (LHS.partialPivLu().solve(RHS)).head(4);
    }
    return M;
}
} // namespace

namespace QI {

std::array<const std::string, 7> TwoPoolModel::varying_names{
    {"PD"s, "T1_m"s, "T2_m"s, "T1_ie"s, "T2_ie"s, "tau_m"s, "f_m"s}};
std::array<const std::string, 2> TwoPoolModel::fixed_names{{"f0"s, "B1"s}};
const QI_ARRAYN(double, 2) TwoPoolModel::fixed_defaults{0.0, 1.0};

TwoPoolModel::TwoPoolModel(const SPGRSequence &s1, const SSFPSequence &s2, const bool scale) :
    spgr{s1}, ssfp{s2}, scale_to_mean{scale} {
    bounds_lo << 1.0, 0.300, 0.010, 0.9, 0.040, 0.025, 0.001;
    bounds_hi << 1.0, 0.800, 0.030, 1.5, 0.150, 0.600, 0.35;
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

bool TwoPoolModel::valid(const QI_ARRAYN(double, NV) & params) const {
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

std::vector<Eigen::ArrayXd> TwoPoolModel::signals(const Eigen::ArrayXd &v,
                                                  const QI_ARRAYN(double, NF) & f) const {
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

Eigen::ArrayXd TwoPoolModel::spgr_signal(const Eigen::ArrayXd &varying,
                                         const QI_ARRAYN(double, NF) & fixed) const {
    const double &  PD    = varying[0];
    const double &  T1_a  = varying[1];
    const double &  T1_b  = varying[3];
    const double &  tau_a = varying[5];
    const double &  f_a   = varying[6];
    const double &  B1    = fixed[1];
    const double &  TR    = spgr.TR;
    Eigen::Matrix2d A, eATR;
    Eigen::Vector2d M0, Mobs;
    Eigen::ArrayXd  signal(spgr.size());
    double          k_ab, k_ba, f_b;
    CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    M0 << f_a, f_b;
    A << ((1. / T1_a) + k_ab), -k_ba, -k_ab, ((1. / T1_b) + k_ba);
    eATR                      = (-TR * A).exp();
    const Eigen::Vector2d RHS = (Eigen::Matrix2d::Identity() - eATR) * M0;
    for (int i = 0; i < spgr.size(); i++) {
        const double a = spgr.FA[i] * B1;
        Mobs = (Eigen::Matrix2d::Identity() - eATR * cos(a)).partialPivLu().solve(RHS * sin(a));
        signal(i) = PD * Mobs.sum();
    }
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}

Eigen::ArrayXd TwoPoolModel::ssfp_signal(const Eigen::ArrayXd &varying,
                                         const QI_ARRAYN(double, NF) & fixed) const {
    Eigen::MatrixXd M       = SSFP2(varying, fixed, ssfp);
    QI_ARRAY(double) signal = M.array().square().colwise().sum().sqrt();
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}

} // End namespace QI
