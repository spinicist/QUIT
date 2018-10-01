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

#include "ThreePoolModel.h"
#include "Helpers.h"

#include <Eigen/Dense>
// There's a bug in matrix log that leads to duplicate symbol definitions if
// this is included in the header file
#include <unsupported/Eigen/MatrixFunctions>

using namespace std::literals;

namespace QI {

std::array<const std::string, ThreePoolModel::NV> ThreePoolModel::varying_names{{"PD"s, "T1_m"s, "T2_m"s, "T1_ie"s, "T2_ie"s, "T1_csf"s, "T2_csf"s, "tau_m"s, "f_m"s, "f_csf"s}};
std::array<const std::string, ThreePoolModel::NF> ThreePoolModel::fixed_names{{"f0"s, "B1"s}};
const QI_ARRAYN(double, ThreePoolModel::NF) ThreePoolModel::fixed_defaults{0.0, 1.0};

/* The inner two-pool model must not be scaled, as scaling will be done once signals are added */
ThreePoolModel::ThreePoolModel(const SequenceType &s, const bool scale) :
    sequence{s}, scale_to_mean{scale}, two_pool{s, false}
{
    bounds_lo << 1.0, 0.300, 0.010, 0.9, 0.040, 3.5, 1.0, 0.025, 0.001, 0.001;
    bounds_hi << 1.0, 0.800, 0.030, 1.5, 0.150, 5.0, 3.5, 0.600, 0.350, 0.999;
}

bool ThreePoolModel::valid(const QI_ARRAYN(double, NV) &params) const {
    // Negative T1/T2 makes no sense
    if ((params[1] <= 0.) || (params[2] <= 0.))
        return false;
    else {
        if ((params[1] < params[3]) && (params[2] < params[4]) && (params[3] < params[5]) &&
            (params[4] < params[6]) && ((params[8] + params[9]) <= 1.0))
            return true;
        else
            return false;
    }
}

auto ThreePoolModel::signal(const Eigen::ArrayXd &varying,
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

auto ThreePoolModel::signal(const Eigen::ArrayXd &varying,
                            const QI_ARRAYN(double, NF) &fixed) const -> Eigen::ArrayXd
{
    Eigen::ArrayXd signals(sequence.size());
    signals.setZero();
    int index = 0;
    for (size_t i = 0; i < sequence.count(); i++) {
        signals.segment(index, sequence.at(i)->size()) = signal(varying, fixed, sequence.at(i));
        index += sequence.at(i)->size();
    }
    return signals;
}

auto ThreePoolModel::signals(const Eigen::ArrayXd &varying,
                             const QI_ARRAYN(double, NF) &fixed) const -> std::vector<Eigen::ArrayXd>
{
    std::vector<Eigen::ArrayXd> signals(sequence.count());
    for (size_t i = 0; i < sequence.count(); i++) {
        signals[i] = signal(varying, fixed, sequence.at(i));
    }
    return signals;
}

Eigen::ArrayXd ThreePoolModel::SSFP1(const double &PD, const double &T1, const double &T2,
                                     const double &f0, const double &B1, const QI::SSFPSequence *s) const
{
    const double E1 = exp(-s->TR/T1);
    const double E2 = exp(-s->TR/T2);
    const Eigen::ArrayXd d = (1 - E1*cos(B1 * s->FA) - (E2*E2)*(E1 - cos(B1 * s->FA)));
    const Eigen::ArrayXd G = sin(B1 * s->FA)*(1 - E1)/d;
    const double a = E2;
    const Eigen::ArrayXd b = E2*(1 - E1)*(1 + cos(B1 * s->FA)) / d;
    const double theta0 = 2.0*M_PI*f0;
    const Eigen::ArrayXd theta = theta0 - s->PhaseInc;
    const double psi = theta0/2.0;
    const Eigen::ArrayXd cos_th = cos(theta);
    const Eigen::ArrayXd sin_th = sin(theta);
    const double cos_psi = cos(psi);
    const double sin_psi = sin(psi);
    const Eigen::ArrayXd re_m = (cos_psi - a*(cos_th*cos_psi - sin_th*sin_psi)) * G / (1.0 - b*cos_th);
    const Eigen::ArrayXd im_m = (sin_psi - a*(cos_th*sin_psi + sin_th*cos_psi)) * G / (1.0 - b*cos_th);
    const Eigen::ArrayXd result = PD * sqrt(re_m.square() + im_m.square());
        // std::cout << PD << " " << T1 << " " << T2 << " " << f0 << " " << B1 << ":" << result.transpose() << std::endl;
    return result;
}

Eigen::ArrayXd ThreePoolModel::spgr_signal(const Eigen::ArrayXd &v,
                                         const QI_ARRAYN(double, NF) &fixed,
                                         const QI::SPGRSequence *s) const
{
    double f_ab = 1. - v[9];
    QI_ARRAYN(double, TwoPoolModel::NV) two_pool_varying;
    two_pool_varying << v[0] * f_ab, v[1], v[2], v[3], v[4], v[7], v[8] / f_ab;
    Eigen::VectorXd m_ab = two_pool.spgr_signal(two_pool_varying, fixed, s);
    Eigen::VectorXd m_c  = SPGRSignal(v[0]*v[9], v[5], fixed[1], *s);
    Eigen::VectorXd signal = m_ab + m_c;
    // std::cout << "SPGR\n" << m_ab.transpose() << "\n" << m_c.transpose() << "\n" << signal.transpose() << std::endl;
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}

Eigen::ArrayXd ThreePoolModel::ssfp_signal(const Eigen::ArrayXd &v,
                                         const QI_ARRAYN(double, NF) &fixed,
                                         const QI::SSFPSequence *s) const
{
    double f_ab = 1. - v[9];
    QI_ARRAYN(double, TwoPoolModel::NV) two_pool_varying;
    two_pool_varying << v[0] * f_ab, v[1], v[2], v[3], v[4], v[7], v[8] / f_ab;
    Eigen::VectorXd m_ab = two_pool.ssfp_signal(two_pool_varying, fixed, s);
    Eigen::VectorXd m_c  = SSFP1(v[0]*v[9], v[5], v[6], fixed[0], fixed[1], s);
    Eigen::VectorXd signal = m_ab + m_c;
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}

Eigen::ArrayXd ThreePoolModel::spgr_signal(const Eigen::ArrayXd &v,
                                             const QI_ARRAYN(double, NF) &fixed,
                                             const QI::SPGREchoSequence *s) const 
{
    double f_ab = 1. - v[9];
    QI_ARRAYN(double, TwoPoolModel::NV) two_pool_varying;
    two_pool_varying << v[0] * f_ab, v[1], v[2], v[3], v[4], v[7], v[8] / f_ab;
    Eigen::VectorXd m_ab = two_pool.spgr_signal(two_pool_varying, fixed, s);
    Eigen::VectorXd m_c  = SPGREchoSignal(v[0]*v[9], v[5], v[6], fixed[1], s);
    Eigen::VectorXd signal = m_ab + m_c;
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}

Eigen::ArrayXd ThreePoolModel::ssfp_signal(const Eigen::ArrayXd &v,
                                         const QI_ARRAYN(double, NF) &fixed,
                                         const QI::SSFPEchoSequence *s) const
{
    double f_ab = 1. - v[9];
    QI_ARRAYN(double, TwoPoolModel::NV) two_pool_varying;
    two_pool_varying << v[0] * f_ab, v[1], v[2], v[3], v[4], v[7], v[8] / f_ab;
    Eigen::VectorXd m_ab = two_pool.ssfp_signal(two_pool_varying, fixed, s);
    Eigen::VectorXd m_c  = SSFP1(v[0]*v[9], v[5], v[6], fixed[0], fixed[1], s) * sqrt(exp(-s->TR / v[6]));
    Eigen::VectorXd signal = m_ab + m_c;
    if (scale_to_mean) {
        signal /= signal.mean();
    }
    return signal;
}


} // End namespace QI
