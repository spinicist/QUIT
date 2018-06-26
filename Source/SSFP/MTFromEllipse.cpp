/*
 *  MTFromEllipseFilter.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "ceres/ceres.h"
#include "MTFromEllipse.h"

namespace QI {

struct EMTCost {
public:
    const Eigen::ArrayXd &G;
    const Eigen::ArrayXd &b;
    const Eigen::ArrayXd flip;
    const Eigen::ArrayXd int_omega2;
    const Eigen::ArrayXd &TR;
    const Eigen::ArrayXd &Trf;
    const double T2_b;
    const double T2_f;
    const double f0_Hz;
    const bool debug;

    template<typename T>
    double debug_print(const T& v) const {
        return v.a;
    }

    template<typename T>
    bool operator() (const T *const p, T* resids) const {
        typedef Eigen::Array<T, Eigen::Dynamic, 1> ArrayXT;
        const T &M0  = p[0];
        const T &F   = p[1];
        const T &k_bf  = p[2];
        const T &T1_f = p[3];
        const T &T1_b = T1_f;

        const ArrayXT E1f = (-TR/T1_f).exp();
        const Eigen::ArrayXd E2_f = (-TR/T2_f).exp();
        const Eigen::ArrayXd E2_fe = (-TR/(2.0*T2_f)).exp();
        const T k_fb = (F > 0.0) ? (k_bf / F) : T(0.0);
        //const double T1_b = 1.0; // Fixed for now
        const ArrayXT E1_b = (-TR/T1_b).exp();
        const ArrayXT fk = (-TR*(k_bf + k_fb)).exp();

        const double G_gauss = (T2_b / sqrt(2.*M_PI))*exp(-pow(2.*M_PI*f0_Hz*T2_b,2) / 2.0);
        const Eigen::ArrayXd WT = M_PI * int_omega2 * G_gauss; // # Product of W and Trf to save a division and multiplication
        const Eigen::ArrayXd fw = (-WT).exp();
        const ArrayXT A = 1.0 + F - fw*E1_b*(F+fk);
        const ArrayXT B = 1.0 + fk*(F-fw*E1_b*(F+1.0));
        const ArrayXT C = F*(1.0-E1_b)*(1.0-fk);

        const ArrayXT denom = (A - B*E1f*cos(flip) - (E2_f*E2_f)*(B*E1f-A*cos(flip)));
        const ArrayXT Gp = M0*E2_fe*(sin(flip)*((1.0-E1f)*B+C))/denom;
        const ArrayXT bp = (E2_f*(A-B*E1f)*(1.0+cos(flip)))/denom;

        Eigen::Map<ArrayXT> r(resids, G.size() + b.size());
        r.head(G.size()) = (G - Gp);
        r.tail(b.size()) = (b - bp);
        if (debug) {
            std::cerr << "M0=" << debug_print(M0) << "\tF=" << debug_print(F)
                      << "\tk_bf=" << debug_print(k_bf) << "\tT1_f=" << debug_print(T1_f) << std::endl;
        }
        return true;
    }
};
template<>
double EMTCost::debug_print(const double &v) const {
    return v;
}

MTFromEllipse::MTFromEllipse(const QI::SSFPMTSequence &s, const double T2, const bool d) :
    m_seq(s), T2_b(T2), debug(d)
{
}

std::vector<float> MTFromEllipse::defaultConsts() const {
    std::vector<float> def(2); // B1, f0
    def[0] = 1.0; def[1] = 0.0;
    return def;
}

MTFromEllipse::TStatus MTFromEllipse::apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                             const TIndex &, // Unused
                             std::vector<TOutput> &outputs, TConst &residual,
                             TInput &resids, TIterations & /* Unused */) const
{
    const double B1 = consts[0];
    const double f0_Hz = consts[1];

    Eigen::Map<const Eigen::ArrayXf> in_G(inputs[0].GetDataPointer(), inputs[0].Size());
    Eigen::Map<const Eigen::ArrayXf> in_a(inputs[1].GetDataPointer(), inputs[1].Size());
    Eigen::Map<const Eigen::ArrayXf> in_b(inputs[2].GetDataPointer(), inputs[2].Size());

    const double scale = in_G.mean();
    Eigen::ArrayXd G = in_G.cast<double>() / scale;
    Eigen::ArrayXd a = in_a.cast<double>();
    Eigen::ArrayXd b = in_b.cast<double>();
    Eigen::ArrayXd T2_fs = (-m_seq.TR / a.log());
    const double T2_f = T2_fs.mean(); // Different TRs so have to average afterwards

    auto *cost = new ceres::AutoDiffCostFunction<EMTCost, ceres::DYNAMIC, 4>(new EMTCost{G, b, m_seq.FA*B1, m_seq.intB1*B1*B1, m_seq.TR, m_seq.Trf, T2_b, T2_f, f0_Hz, debug}, G.size() + b.size());
    ceres::LossFunction *loss = new ceres::HuberLoss(1.0);
    Eigen::Array<double, 4, 1> p; p << 15.0, 0.1, 2.5, 1.0;
    ceres::Problem problem;
    problem.AddResidualBlock(cost, loss, p.data());
    problem.SetParameterLowerBound(p.data(), 0, 0.1);
    problem.SetParameterUpperBound(p.data(), 0, 20.0);
    problem.SetParameterLowerBound(p.data(), 1, 1e-6);
    problem.SetParameterUpperBound(p.data(), 1, 0.2 - 1e-6);
    problem.SetParameterLowerBound(p.data(), 2, 0.1);
    problem.SetParameterUpperBound(p.data(), 2, 20.0);
    problem.SetParameterLowerBound(p.data(), 3, 0.05);
    problem.SetParameterUpperBound(p.data(), 3, 5.0);
    ceres::Solver::Options options;
    ceres::Solver::Summary summary;
    options.max_num_iterations = 100;
    options.function_tolerance = 1e-7;
    options.gradient_tolerance = 1e-8;
    options.parameter_tolerance = 1e-3;
    if (!debug) options.logging_type = ceres::SILENT;
    ceres::Solve(options, &problem, &summary);
    if (!summary.IsSolutionUsable()) {
        return std::make_tuple(false, summary.FullReport());
    } else if (debug) {
        std::cout << summary.FullReport() << std::endl;
    }
    outputs[0] = p[0] * scale;
    outputs[1] = p[1] / (1 + p[1]); // f = F/(1+F);
    outputs[2] = p[2];
    outputs[3] = p[3];
    outputs[4] = T2_f;
    residual = summary.final_cost;
    if (resids.Size() > 0) {
        assert(resids.Size() == (G.size() + a.size() + b.size()));
        std::vector<double> r_temp(G.size() + a.size() + b.size());
        problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
        for (int i = 0; i < G.size(); i++)
            resids[i] = r_temp[i];
        Eigen::ArrayXd as = (-m_seq.TR / T2_f).exp();
        for (int i = 0; i < a.size(); i++) {
            resids[i + G.size()] = as[i] - a[i];
        }
        for (int i = 0; i < b.size(); i++) {
            resids[i + G.size() + a.size()] = r_temp[i + G.size()];
        }
    }
    return std::make_tuple(true, "");
}

} // End namespace QI