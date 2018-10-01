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

#include "TwoPoolModel.h"
#include "Helpers.h"

// There's a bug in matrix log that leads to duplicate symbol definitions if
// this is included in the header file
#include <unsupported/Eigen/MatrixFunctions>

using namespace std::literals;

namespace QI {

std::array<const std::string, 7> TwoPoolModel::varying_names{{"PD"s, "T1_m"s, "T2_m"s, "T1_ie"s, "T2_ie"s, "tau_m"s, "f_m"s}};
std::array<const std::string, 2> TwoPoolModel::fixed_names{{"f0"s, "B1"s}};
const QI_ARRAYN(double, 2) TwoPoolModel::fixed_defaults{0.0, 1.0};

TwoPoolModel::TwoPoolModel(const SequenceType &s, const bool scale) :
    sequence{s}, scale_to_mean{scale}
{
    bounds_lo << 1.0, 0.300, 0.010, 0.9, 0.040, 0.025, 0.001;
    bounds_hi << 1.0, 0.800, 0.030, 1.5, 0.150, 0.600, 0.35;
}

bool TwoPoolModel::valid(const QI_ARRAYN(double, NV) &params) const {
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

auto TwoPoolModel::signal(const Eigen::ArrayXd &varying,
                          const QI_ARRAYN(double, NF) &fixed,
                          const QI::SequenceBase *s) const -> Eigen::ArrayXd
{
    const auto spgr = dynamic_cast<const QI::SPGRSequence *>(s);
    if (spgr) return spgr_signal(varying, fixed, spgr);
    const auto ssfp = dynamic_cast<const QI::SSFPSequence *>(s);
    if (ssfp) return ssfp_signal(varying, fixed, ssfp);
    const auto spgre = dynamic_cast<const QI::SPGREchoSequence *>(s);
    if (spgre) return spgr_signal(varying, fixed, spgre);
    const auto ssfpe = dynamic_cast<const QI::SSFPEchoSequence *>(s);
    if (ssfpe) return ssfp_signal(varying, fixed, ssfpe);
    QI_FAIL("Given pointer was not to SPGR/SPGREcho/SSFP/SSFPEcho");
}

auto TwoPoolModel::signal(const Eigen::ArrayXd &varying,
                          const QI_ARRAYN(double, NF) &fixed) const -> Eigen::ArrayXd
{
    Eigen::ArrayXd signals(sequence.size());
    int index = 0;
    for (size_t i = 0; i < sequence.count(); i++) {
        signals.segment(index, sequence.at(i)->size()) = signal(varying, fixed, sequence.at(i));
        index += sequence.at(i)->size();
    }
    return signals;
}

auto TwoPoolModel::signals(const Eigen::ArrayXd &varying,
                           const QI_ARRAYN(double, NF) &fixed) const -> std::vector<Eigen::ArrayXd>
{
    std::vector<Eigen::ArrayXd> signals(sequence.count());
    for (size_t i = 0; i < sequence.count(); i++) {
        signals[i] = signal(varying, fixed, sequence.at(i));
    }
    return signals;
}

Eigen::ArrayXd TwoPoolModel::spgr_signal(const Eigen::ArrayXd &varying,
                                         const QI_ARRAYN(double, NF) &fixed,
                                         const QI::SPGRSequence *s) const
{
    const double &PD    = varying[0];
    const double &T1_a  = varying[1];
    const double &T1_b  = varying[3];
    const double &tau_a = varying[5];
    const double &f_a   = varying[6];
    const double &B1 = fixed[1];
    const double &TR = s->TR;
    Eigen::Matrix2d A, eATR;
    Eigen::Vector2d M0, Mobs;
    Eigen::ArrayXd signal(s->size());
    double k_ab, k_ba, f_b;
    CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    M0 << f_a, f_b;
    A << ((1./T1_a) + k_ab),            -k_ba,
                     -k_ab, ((1./T1_b) + k_ba);
    eATR = (-TR*A).exp();
    const Eigen::Vector2d RHS = (Eigen::Matrix2d::Identity() - eATR) * M0;
    for (int i = 0; i < s->size(); i++) {
        const double a = s->FA[i] * B1;
        Mobs = (Eigen::Matrix2d::Identity() - eATR*cos(a)).partialPivLu().solve(RHS * sin(a));
        signal(i) = PD * Mobs.sum();
    }
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;

}

Eigen::MatrixXd TwoPoolModel::SSFP2(const Eigen::ArrayXd &varying,
                                    const QI_ARRAYN(double, NF) &fixed,
                                    const QI::SSFPSequence *s) const
{
    const double &PD    = varying[0];
    const double &T1_a  = varying[1];
    const double &T2_a  = varying[2];
    const double &T1_b  = varying[3];
    const double &T2_b  = varying[4];
    const double &tau_a = varying[5];
    const double &f_a   = varying[6];
    const double &f0 = fixed[0];
    const double &B1 = fixed[1];
    const double &TR = s->TR;
    const double E1_a = exp(-TR/T1_a);
    const double E1_b = exp(-TR/T1_b);
    const double E2_a = exp(-TR/T2_a);
    const double E2_b = exp(-TR/T2_b);
    double f_b, k_ab, k_ba;
    CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    const double E_ab = exp(-TR*k_ab/f_b);
    const double K1 = E_ab*f_b+f_a;
    const double K2 = E_ab*f_a+f_b;
    const double K3 = f_a*(1-E_ab);
    const double K4 = f_b*(1-E_ab);
    const Eigen::ArrayXd alpha = B1 * s->FA;
    const Eigen::ArrayXd theta = s->PhaseInc + 2.*M_PI*f0*TR;

    Eigen::MatrixXd M(4, s->size());
    Eigen::Matrix6d LHS;
    Eigen::Vector6d RHS;
    RHS << 0, 0, 0, 0, -E1_b*K3*f_b + f_a*(-E1_a*K1 + 1), -E1_a*K4*f_a + f_b*(-E1_b*K2 + 1);

    for (int i = 0; i < s->size(); i++) {
        const double ca = cos(alpha[i]);
        const double sa = sin(alpha[i]);
        const double cta = cos(theta[i]);
        const double ctb = cos(theta[i]);
        const double sta = sin(theta[i]);
        const double stb = sin(theta[i]);

        LHS << -E2_a*K1*cta + ca, -E2_b*K3*cta, E2_a*K1*sta, E2_b*K3*sta, sa, 0,
            -E2_a*K4*ctb, -E2_b*K2*ctb + ca, E2_a*K4*stb, E2_b*K2*stb, 0, sa,
            -E2_a*K1*sta, -E2_b*K3*sta, -E2_a*K1*cta + 1, -E2_b*K3*cta, 0, 0,
            -E2_a*K4*stb, -E2_b*K2*stb, -E2_a*K4*ctb, -E2_b*K2*ctb + 1, 0, 0,
            -sa, 0, 0, 0, -E1_a*K1 + ca, -E1_b*K3,
                0, -sa, 0, 0, -E1_a*K4, -E1_b*K2 + ca;
        M.col(i).noalias() = PD * (LHS.partialPivLu().solve(RHS)).head(4);
    }
    return M;
}

Eigen::ArrayXd TwoPoolModel::ssfp_signal(const Eigen::ArrayXd &varying,
                                         const QI_ARRAYN(double, NF) &fixed,
                                         const QI::SSFPSequence *s) const
{
    Eigen::MatrixXd M = SSFP2(varying, fixed, s);
    QI_ARRAY(double) signal = M.array().square().colwise().sum().sqrt();
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}

Eigen::ArrayXd TwoPoolModel::spgr_signal(const Eigen::ArrayXd &varying,
                                             const QI_ARRAYN(double, NF) &fixed,
                                             const QI::SPGREchoSequence *s) const 
{
    const double &PD    = varying[0];
    const double &T1_a  = varying[1];
    const double &T2_a  = varying[2];
    const double &T1_b  = varying[3];
    const double &T2_b  = varying[4];
    const double &tau_a = varying[5];
    const double &f_a   = varying[6];
    const double &B1 = fixed[1];
    const double &TR = s->TR;
    const double &TE = s->TE;
    Eigen::Matrix2d A, eATR;
    Eigen::Vector2d M0, Mz, Mxy;
    double k_ab, k_ba, f_b;
    CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    M0 << f_a, f_b;
    A << -(1./T1_a) - k_ab,              k_ba,
                      k_ab, -(1./T1_b) - k_ba;
    eATR = (TR*A).exp();
    Eigen::Matrix2d echo = Eigen::Matrix2d::Zero();
    echo(0,0) = exp(-TE/T2_a); // T2' absorbed into PD as it effects both components equally
    echo(1,1) = exp(-TE/T2_b);
    const Eigen::Vector2d RHS = (Eigen::Matrix2d::Identity() - eATR) * M0;
    Eigen::VectorXd signal(s->size());
    for (int i = 0; i < s->size(); i++) {
        const double a = s->FA[i] * B1;
        Mz.noalias() = (Eigen::Matrix2d::Identity() - eATR*cos(a)).partialPivLu().solve(RHS);
        Mxy.noalias() = echo * Mz * sin(a);
        signal(i) = PD * (Mxy(0) + Mxy(1));
    }
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}

Eigen::ArrayXd TwoPoolModel::ssfp_signal(const Eigen::ArrayXd &varying,
                                         const QI_ARRAYN(double, NF) &fixed,
                                         const QI::SSFPEchoSequence *s) const
{
    Eigen::MatrixXd M = SSFP2(varying, fixed, s);
    const double &T2_a  = varying[2];
    const double &T2_b  = varying[4];
    const double &tau_a = varying[5];
    const double &f_a   = varying[6];
    const double &f0    = fixed[0];
    const double &TR = s->TR;
    const double sE2_a = exp(-TR/(2.*T2_a));
    const double sE2_b = exp(-TR/(2.*T2_b));
    double f_b, k_ab, k_ba;
    CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
    const double sqrtE_ab = exp(-TR*k_ab/(2*f_b));
    const double K1 = sqrtE_ab*f_b+f_a;
    const double K2 = sqrtE_ab*f_a+f_b;
    const double K3 = f_a*(1-sqrtE_ab);
    const double K4 = f_b*(1-sqrtE_ab);
    const double cp = cos(M_PI*f0*TR);
    const double sp = sin(M_PI*f0*TR);

    Eigen::Matrix<double, 4, 4> echo;
    echo << sE2_a*K1*cp, sE2_b*K3*cp, -sE2_a*K1*sp, -sE2_b*K3*sp,
            sE2_a*K4*cp, sE2_b*K2*cp, -sE2_a*K4*sp, -sE2_b*K2*sp,
            sE2_a*K1*sp, sE2_b*K3*sp,  sE2_a*K1*cp,  sE2_b*K3*cp,
            sE2_a*K4*sp, sE2_b*K2*sp,  sE2_a*K4*cp,  sE2_b*K2*cp;

    const Eigen::MatrixXd Me = echo*M;
    QI_ARRAY(double) signal = Me.array().square().colwise().sum().sqrt();
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}


} // End namespace QI
