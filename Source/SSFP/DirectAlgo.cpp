/*
 *  DirectAlgo.cpp
 *
 *  Copyright (c) 2016, 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "DirectAlgo.h"
#include "EllipseHelpers.h"
#include "Fit.h"
#include "ceres/ceres.h"

namespace QI {

struct DirectCost {
public:
    const Eigen::ArrayXcd &data;
    const double TR;
    const Eigen::ArrayXd &phi;
    const bool debug;

    template<typename T>
    bool operator() (T const* const* p, T* resids) const {
        typedef Eigen::Array<T, Eigen::Dynamic, 1> ArrayXT;

        const T &G = p[0][0];
        const T &a = p[0][1];
        const T &b = p[0][2];
        const T &th0 = p[0][3];
        const T &psi0 = p[0][4];

        if (b < 2.*a/(1. + a*a)) {
            ArrayXT m = EllipseToSignal(G, a, b, th0, psi0, TR, phi);
            Eigen::Map<ArrayXT> r(resids, data.size()*2);
            r.head(data.size()) = m.head(data.size()) - data.real();
            r.tail(data.size()) = m.tail(data.size()) - data.imag();
            if (debug) {
                std::cout << "*** ELLIPSE COST *** "
                    << "G " << G << " a " << a << " b " << b << " th0 " << th0 << " psi0 " << psi0 << std::endl;
            }
            return true;
        } else {
            if (debug) {
                std::cout << "*** ELLIPSE COST NOT EVALUATED - ELLIPSE WAS HORIZONTAL *** " 
                    << "G " << G << " a " << a << " b " << b << " th0 " << th0 << " psi0 " << psi0 << std::endl;
            }
            return false;
        }
    }
};

Eigen::ArrayXd DirectAlgo::apply_internal(const Eigen::ArrayXcf &indata, const double flip, const double TR, const Eigen::ArrayXd &phi, const bool debug, float &residual) const {
    Eigen::ArrayXcd data = indata.cast<std::complex<double>>();
    const double scale = data.abs().maxCoeff();
    data /= scale;

    std::complex<double> c_mean = data.mean();

    auto *cost = new ceres::DynamicAutoDiffCostFunction<DirectCost>(new DirectCost{data, TR, phi, debug});
    cost->AddParameterBlock(5);
    cost->SetNumResiduals(data.size()*2);
    ceres::LossFunction *loss = new ceres::HuberLoss(1.0);
    Eigen::Array<double, 5, 1> p, best_p;
    ceres::Problem problem;
    problem.AddResidualBlock(cost, loss, p.data());
    const double not_zero = 1.0e-6;
    const double not_one  = 1.0 - not_zero;
    const double max_a = exp(-TR / 4.3); // Set a sensible maximum on T2
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
    options.parameter_tolerance = 1e-4;
    if (!debug) options.logging_type = ceres::SILENT;

    // Calculate a sensible guess for a/b using T1/T2 of grey matter
    Eigen::Array3d Gab = EllipseGab(1000., 50., TR, flip);
    double best_resid = std::numeric_limits<double>::infinity();
    for (auto &th0_start : {0., -M_PI, M_PI}) {
        const double psi0_est   = arg(c_mean / std::polar(1.0, th0_start/2));
        p << abs(c_mean), Gab[1], Gab[2], th0_start, psi0_est;
        ceres::Solve(options, &problem, &summary);
        if (debug || !summary.IsSolutionUsable()) {
            std::cout << summary.FullReport() << std::endl;
        }
        if (summary.final_cost < best_resid) {
            if (debug) std::cout << "*** Final cost is an improvement. Saving" << std::endl;
            best_p = p;
            best_resid = summary.final_cost;
        }
    }
    residual = best_resid;
    best_p[0] *= scale;
    best_p[3] = std::fmod(best_p[3] + 3*M_PI, 2*M_PI) - M_PI;
    best_p[4] = std::fmod(best_p[4] + 3*M_PI, 2*M_PI) - M_PI;
    if (debug) {
        std::cout << "Finished Direct Algo, p: " << p.transpose() << std::endl;
    }
    return best_p;
};

} // End namespace QI