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
    bool operator() (const T *const p, T* resids) const {
        typedef Eigen::Array<T, Eigen::Dynamic, 1> ArrayXT;

        const T &G = p[0];
        const T &a = p[1];
        const T &b = p[2];
        const T &th0 = p[3];
        const T &psi0 = p[4];

        if (b < 2.*a/(1. + a*a)) {
            ArrayXT m = EllipseToSignal(G, a, b, th0, psi0, TR, phi);
            Eigen::Map<ArrayXT> r(resids, data.size()*2);
            r.head(data.size()) = m.head(data.size()) - data.real();
            r.tail(data.size()) = m.tail(data.size()) - data.imag();
            Eigen::IOFormat numpy(4, 0, ", ", ",", "[", "]", "np.array(", ")");
            if (debug) {
                Eigen::ArrayXd md(2*data.size()); for (int i = 0; i < 2*data.size(); i++) { md[i] = m[i].a; }
                std::cout 
                    << "{ 'P': { 'G': " << G.a << ", 'a': " << a.a << " , 'b': " << b.a << " , 'th0': " << th0.a << " , 'psi0': " << psi0.a << "},"
                    << " 'M': { 'real': " << md.head(data.size()).transpose().format(numpy) << ", 'imag': " << md.tail(data.size()).transpose().format(numpy) << "},"
                    << " 'D': { 'real': " << data.real().transpose().format(numpy) << ", 'imag': " << data.imag().transpose().format(numpy) << "}}," << std::endl;
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

    template<>
    bool DirectCost::operator() <double> (const double *const p, double *resids) const {
        const double &G = p[0];
        const double &a = p[1];
        const double &b = p[2];
        const double &th0 = p[3];
        const double &psi0 = p[4];
        
        if (b < 2.*a/(1. + a*a)) {
            Eigen::ArrayXd m = EllipseToSignal(G, a, b, th0, psi0, TR, phi);
            Eigen::Map<Eigen::ArrayXd> r(resids, data.size()*2);
            r.head(data.size()) = m.head(data.size()) - data.real();
            r.tail(data.size()) = m.tail(data.size()) - data.imag();
            if (debug) {

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


Eigen::ArrayXd DirectAlgo::apply_internal(const Eigen::ArrayXcf &indata,
                                          const double flip, const double TR, const Eigen::ArrayXd &phi,
                                          const bool debug, float &residual) const {
    Eigen::ArrayXcd data = indata.cast<std::complex<double>>();
    const double scale = data.abs().maxCoeff();
    data /= scale;
    std::complex<double> c_mean = data.mean();

    auto *cost = new ceres::AutoDiffCostFunction<DirectCost, ceres::DYNAMIC, 5>(new DirectCost{data, TR, phi, debug}, data.size()*2);
    ceres::LossFunction *loss = new ceres::HuberLoss(1.0);
    Eigen::Array<double, 5, 1> p;
    ceres::Problem problem;
    problem.AddResidualBlock(cost, loss, p.data());
    const double not_zero = 1.0e-6;
    const double not_one  = 1.0 - not_zero;
    const double max_a = exp(-TR / 4.3); // Set a sensible maximum on T2
    if (debug) std::cout << "max_a : " << max_a << std::endl;
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

    // Calculate a sensible guess for a/b using T1/T2 of grey matter
    Eigen::Array3d Gab = EllipseGab(1.0, 0.05, TR, flip);
    std::array<double, 3> costs, start_th0{-M_PI, 0., M_PI};
    double th0, psi0, best_cost = std::numeric_limits<double>::infinity();
    for (const auto &th0_try : {-M_PI, 0., M_PI}) {
        const double psi0_try = arg(c_mean / std::polar(1.0, th0_try/2));
        p << abs(c_mean), Gab[1], Gab[2], th0_try, psi0_try;
        double cost;
        problem.Evaluate(ceres::Problem::EvaluateOptions(), &cost, NULL, NULL, NULL);
        if (debug) std::cout << "th0 = " << th0_try << ", cost was " << cost;
        if (cost < best_cost) {
            best_cost = cost;
            th0 = th0_try;
            psi0 = psi0_try;
            if (debug) std::cout << ", saving for start" << std::endl;
        } else {
            if (debug) std::cout << std::endl;
        }
    }
    p << abs(c_mean), Gab[1], Gab[2], th0, psi0;
    if (debug) std::cout << "Starting p: " << p.transpose() << std::endl;
    ceres::Solve(options, &problem, &summary);
    if (debug || !summary.IsSolutionUsable()) {
        std::cout << summary.FullReport() << std::endl;
    }
    residual = summary.final_cost;
    p[0] *= scale;
    p[3] = std::fmod(p[3] + 3*M_PI, 2*M_PI) - M_PI;
    p[4] = std::fmod(p[4] + 3*M_PI, 2*M_PI) - M_PI;
    return p;
};

} // End namespace QI