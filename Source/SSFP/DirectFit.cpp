/*
 *  DirectFit.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "DirectFit.h"
#include "Macro.h"
#include "ceres/ceres.h"

namespace QI {

struct DirectCost : ModelCost<EllipseModel> {
    using ModelCost::ModelCost; // Inherit constructor

    template<typename T>
    bool operator() (const T *const vin, T* rin) const {
        const Eigen::Map<const QI_ARRAYN(T, EllipseModel::NV)> v(vin);

        // Check if ellipse has gone horizontal
        const T &a = v[1];
        const T &b = v[2];
        if (b < 2.*a/(1. + a*a)) {
            QI_ARRAY(T) m = model.signal(v, fixed, sequence);
            Eigen::Map<QI_ARRAY(T)> r(rin, data.size()*2);
            r.head(data.size()) = m.head(data.size()) - data.real();
            r.tail(data.size()) = m.tail(data.size()) - data.imag();
            // if (debug) {
            //     Eigen::IOFormat numpy(4, 0, ", ", ",", "[", "]", "np.array(", ")");
            //     Eigen::ArrayXd md(2*data.size()); for (int i = 0; i < 2*data.size(); i++) { md[i] = m[i].a; }
            //     std::cout 
            //         << "{ 'P': { 'G': " << G.a << ", 'a': " << a.a << " , 'b': " << b.a << " , 'th0': " << th0.a << " , 'psi0': " << psi0.a << "},"
            //         << " 'M': { 'real': " << md.head(data.size()).transpose().format(numpy) << ", 'imag': " << md.tail(data.size()).transpose().format(numpy) << "},"
            //         << " 'D': { 'real': " << data.real().transpose().format(numpy) << ", 'imag': " << data.imag().transpose().format(numpy) << "}}," << std::endl;
            // }
            return true;
        } else {
            // if (debug) {
            //     std::cout << "*** ELLIPSE COST NOT EVALUATED - ELLIPSE WAS HORIZONTAL *** " 
            //         << "G " << G << " a " << a << " b " << b << " th0 " << th0 << " psi0 " << psi0 << std::endl;
            // }
            return false;
        }
    }
};

QI::FitReturnType DirectFit::fit(const std::vector<Eigen::ArrayXcd> &inputs,
                                 const Eigen::ArrayXd &fixed, QI_ARRAYN(OutputType, EllipseModel::NV) &p,
                                 ResidualType &residual, std::vector<Eigen::ArrayXcd> &/*Unused*/, FlagType &iterations,
                                 const int /* Unused */) const
{
    const double scale = inputs[0].abs().maxCoeff();
    const Eigen::ArrayXcd data = inputs[0] / scale;
    const std::complex<double> c_mean = data.mean();

    using AutoCost = ceres::AutoDiffCostFunction<DirectCost, ceres::DYNAMIC, EllipseModel::NV>;
    auto *cost = new DirectCost(model, sequence, fixed, data);
    auto *auto_cost = new AutoCost(cost, data.rows() * 2);
    ceres::LossFunction *loss = new ceres::HuberLoss(1.0);
    ceres::Problem problem;
    problem.AddResidualBlock(auto_cost, loss, p.data());
    const double not_zero = 1.0e-6;
    const double not_one  = 1.0 - not_zero;
    const double max_a = exp(-sequence->TR / 5.0); // Set a sensible maximum on T2
    problem.SetParameterLowerBound(p.data(), 0, not_zero); problem.SetParameterUpperBound(p.data(), 0, not_one);
    problem.SetParameterLowerBound(p.data(), 1, not_zero); problem.SetParameterUpperBound(p.data(), 1, max_a);
    problem.SetParameterLowerBound(p.data(), 2, not_zero); problem.SetParameterUpperBound(p.data(), 2, not_one);
    problem.SetParameterLowerBound(p.data(), 3, -2.*M_PI); problem.SetParameterUpperBound(p.data(), 3, 2.*M_PI);
    problem.SetParameterLowerBound(p.data(), 4, -2.*M_PI); problem.SetParameterUpperBound(p.data(), 4, 2.*M_PI);
    ceres::Solver::Options options;
    ceres::Solver::Summary summary;
    options.max_num_iterations = 50;
    options.function_tolerance = 1e-5;
    options.gradient_tolerance = 1e-6;
    options.parameter_tolerance = 1e-3;
    options.logging_type = ceres::SILENT;
    double th0, psi0, best_cost = std::numeric_limits<double>::infinity();
    for (const auto &th0_try : {-M_PI, 0., M_PI}) {
        const double psi0_try = arg(c_mean / std::polar(1.0, th0_try/2));
        p << abs(c_mean), 0.5, 0.5, th0_try, psi0_try;
        double cost;
        problem.Evaluate(ceres::Problem::EvaluateOptions(), &cost, NULL, NULL, NULL);
        if (cost < best_cost) {
            best_cost = cost;
            th0 = th0_try;
            psi0 = psi0_try;
            iterations = summary.iterations.size();
            residual = summary.final_cost * scale;
        }
    }
    p << abs(c_mean), 0.5, 0.5, th0, psi0;
    ceres::Solve(options, &problem, &summary);
    if (!summary.IsSolutionUsable()) {
        return std::make_tuple(false, summary.FullReport());
    }
    residual = summary.final_cost;
    p[0] *= scale;
    p[3] = std::fmod(p[3] + 3*M_PI, 2*M_PI) - M_PI;
    p[4] = std::fmod(p[4] + 3*M_PI, 2*M_PI) - M_PI;
    return std::make_tuple(true, "");
}

} // End namespace QI