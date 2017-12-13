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
    const Eigen::ArrayXd &flip;
    const Eigen::ArrayXd &int_omega2;
    const Eigen::ArrayXd &TR;
    const Eigen::ArrayXd &Trf;
    const double T2r;
    const double T2f;
    const double f0_Hz;
    const bool debug;

    template<typename T>
    bool operator() (T const* const* p, T* resids) const {
        typedef Eigen::Array<T, Eigen::Dynamic, 1> ArrayXT;

        const T &M0  = p[0][0];
        const T &F   = p[0][1];
        const T &kf  = p[0][2];
        const T &T1f = p[0][3];
        const T &T1r = T1f; //p[0][3];
        //const double &T2r = 12.0e-6; //p[0][4];

        const ArrayXT E1f = (-TR/T1f).exp();
        const Eigen::ArrayXd E2f = (-TR/T2f).exp();
        const Eigen::ArrayXd E2f_echo = (-TR/(2.0*T2f)).exp();
        const T kr = (F > 0.0) ? (kf / F) : T(0.0);
        //const double T1r = 1.0; // Fixed for now
        const ArrayXT E1r = (-TR/T1r).exp();
        const ArrayXT fk = (-TR*(kf + kr)).exp();

        //const double T2r = 25.e-6; // Fixed for now;
        const double G_gauss = (T2r / sqrt(2.*M_PI))*exp(-pow(2.*M_PI*f0_Hz*T2r,2) / 2.0);
        const Eigen::ArrayXd WT = M_PI * int_omega2 * G_gauss; // # Product of W and Trf to save a division and multiplication
        const Eigen::ArrayXd fw = (-WT).exp();
        const ArrayXT A = 1.0 + F - fw*E1r*(F+fk);
        const ArrayXT B = 1.0 + fk*(F-fw*E1r*(F+1.0));
        const ArrayXT C = F*(1.0-E1r)*(1.0-fk);

        const ArrayXT denom = (A - B*E1f*cos(flip) - (E2f*E2f)*(B*E1f-A*cos(flip)));
        const ArrayXT Gp = M0*E2f_echo*(sin(flip)*((1.0-E1f)*B+C))/denom;
        const ArrayXT bp = (E2f*(A-B*E1f)*(1.0+cos(flip)))/denom;

        Eigen::Map<ArrayXT> r(resids, G.size() + b.size());
        r.head(G.size()) = (G - Gp);
        r.tail(b.size()) = (b - bp);
        if (debug) {
            std::cerr << "M0=" << M0 << "\nF=" << F << "\nkf=" << kf << "\nT1f=" << T1f << "\nT2f=" << T2f << "\nf0_Hz=" << f0_Hz << "\n"
                      << "flip:" << flip.transpose() << "\n"
                      << "int: " << int_omega2.transpose() << std::endl;
        }
        return true;
    }
};

MTFromEllipse::MTFromEllipse(const Eigen::ArrayXd &f, const Eigen::ArrayXd &iB, const Eigen::ArrayXd &tr, const Eigen::ArrayXd &trf, const double T2, const bool d) :
    flips(f), intB1(iB), TRs(tr), TRFs(trf), T2r(T2), debug(d)
{
}

std::vector<float> MTFromEllipse::defaultConsts() const {
    std::vector<float> def(2); // B1, f0
    def[0] = 1.0; def[1] = 0.0;
    return def;
}

bool MTFromEllipse::apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                    std::vector<TOutput> &outputs, TConst &residual,
                    TInput &resids, TIters &its) const
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
    Eigen::ArrayXd T2fs = (-TRs.cast<double>() / a.log());
    const double T2f = T2fs.mean(); // Different TRs so have to average afterwards

    auto *cost = new ceres::DynamicAutoDiffCostFunction<EMTCost>(new EMTCost{G, b, flips*B1, intB1*B1*B1, TRs, TRFs, T2r, T2f, f0_Hz, debug});
    cost->AddParameterBlock(4);
    cost->SetNumResiduals(G.size() + b.size());
    ceres::LossFunction *loss = new ceres::HuberLoss(1.0);
    Eigen::Array<double, 6, 1> p; p << 15.0, 0.05, 5.0, 1.0;
    ceres::Problem problem;
    problem.AddResidualBlock(cost, loss, p.data());
    const double not_zero = std::nextafter(0.0, 1.0);
    const double not_one  = std::nextafter(1.0, 0.0);
    problem.SetParameterLowerBound(p.data(), 0, 0.1);
    problem.SetParameterUpperBound(p.data(), 0, 20.0);
    problem.SetParameterLowerBound(p.data(), 1, 1e-6);
    problem.SetParameterUpperBound(p.data(), 1, 0.2 - 1e-6);
    problem.SetParameterLowerBound(p.data(), 2, 0.1);
    problem.SetParameterUpperBound(p.data(), 2, 10.0);
    problem.SetParameterLowerBound(p.data(), 3, 0.5);
    problem.SetParameterUpperBound(p.data(), 3, 5.0);
    ceres::Solver::Options options;
    ceres::Solver::Summary summary;
    options.max_num_iterations = 100;
    options.function_tolerance = 1e-7;
    options.gradient_tolerance = 1e-8;
    options.parameter_tolerance = 1e-6;
    if (!debug) options.logging_type = ceres::SILENT;
    ceres::Solve(options, &problem, &summary);
    if (!summary.IsSolutionUsable()) {
        std::cerr << summary.FullReport() << std::endl;
        std::cerr << "Parameters: " << p.transpose() << " T2f: " << T2f << " B1: " << B1 << std::endl;
        std::cerr << "G: " << G.transpose() << std::endl;
        std::cerr << "a: " << a.transpose() << std::endl;
        std::cerr << "b: " << b.transpose() << std::endl;
        return false;
    } else if (debug) {
        std::cout << summary.FullReport() << std::endl;
    }
    outputs[0] = p[0] * scale;
    outputs[1] = p[1];
    outputs[2] = p[2];
    outputs[3] = p[3];
    outputs[4] = T2f;
    residual = summary.final_cost;
    if (resids.Size() > 0) {
        assert(resids.Size() == data.size());
        std::vector<double> r_temp(G.size() + a.size() + b.size());
        problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
        for (int i = 0; i < G.size(); i++)
            resids[i] = r_temp[i];
        Eigen::ArrayXd as = (-TRs.cast<double>() / T2f).exp();
        for (int i = 0; i < a.size(); i++) {
            resids[i + G.size()] = as[i] - a[i];
        }
        for (int i = 0; i < b.size(); i++) {
            resids[i + G.size() + a.size()] = r_temp[i + G.size()];
        }
    }
    return true;
}

} // End namespace QI