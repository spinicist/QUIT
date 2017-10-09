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

#include "QI/Ellipse/DirectAlgo.h"
#include "QI/Ellipse/EllipseHelpers.h"
#include "QI/Fit.h"
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
        const T &f0 = p[0][3];
        const T &psi0 = p[0][4];

        ArrayXT m = EllipseToSignal(G, a, b, f0, psi0, TR, phi);
        Eigen::Map<ArrayXT> r(resids, data.size()*2);
        r.head(data.size()) = m.head(data.size()) - data.real();
        r.tail(data.size()) = m.tail(data.size()) - data.imag();
        if (debug) {
            std::cout << "*** COST ***\n"
                << "G " << G << " a " << a << " b " << b << " f0 " << f0 << " psi0 " << psi0;
        }
        return true;
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

    // Get as estimate of f0 and psi0
    // Do a linear regression of phi against unwrapped phase diff, intercept is theta0
    Eigen::VectorXd Y = QI::Unwrap((data / c_mean - std::complex<double>(1.0, 0.0)).arg());
    Eigen::MatrixXd X(Y.rows(), 2);
    X.col(0) = phi - M_PI;
    X.col(1).setOnes();
    const Eigen::VectorXd b = QI::LeastSquares(X, Y);
    const double theta0_est = fmod(arg((data[0] / c_mean) - std::complex<double>(1.0, 0.0)), M_PI/2);//b[1];
    const double psi0_est   = arg(c_mean / std::polar(1.0, theta0_est/2));
    const double f0_est     = theta0_est / (2.0 * M_PI * TR);
    if (debug) {
        std::cerr << "X\n" << X.transpose() << "\nY\n" << Y.transpose() << std::endl;
        std::cerr << "b " << b.transpose() << " theta0_est " << theta0_est << " psi0_est " << psi0_est << std::endl;
        std::cerr << "arg(c_mean) " << std::arg(c_mean) << std::endl;
        std::cerr << "old theta0_est " << arg((data[0] / c_mean) - std::complex<double>(1.0, 0.0)) << std::endl;
    }
    
    Eigen::Array<double, 5, 1> p;

    // Calculate a sensible guess for a/b using T1/T2 of grey matter
    Eigen::Array3d Gab = EllipseGab(1000., 50., TR, flip);

    p << abs(c_mean), Gab[1], Gab[2], f0_est, psi0_est;
    // std::cout << "Start p: " << p.transpose() << std::endl;
    ceres::Problem problem;
    problem.AddResidualBlock(cost, NULL, p.data());
    const double not_zero = 1.0e-6;
    const double not_one  = 1.0 - not_zero;
    const double max_a = exp(-TR / 4.3); // Set a sensible maximum on T2
    problem.SetParameterLowerBound(p.data(), 0, not_zero); problem.SetParameterUpperBound(p.data(), 0, not_one);
    problem.SetParameterLowerBound(p.data(), 1, not_zero); problem.SetParameterUpperBound(p.data(), 1, max_a);
    problem.SetParameterLowerBound(p.data(), 2, not_zero); problem.SetParameterUpperBound(p.data(), 2, not_one);
    problem.SetParameterLowerBound(p.data(), 3, -0.5/TR + 1.0e-6); problem.SetParameterUpperBound(p.data(), 3,  0.5/TR - 1.0e-6);
    problem.SetParameterLowerBound(p.data(), 4, -M_PI); problem.SetParameterUpperBound(p.data(), 4,  M_PI);
    ceres::Solver::Options options;
    ceres::Solver::Summary summary;
    options.max_num_iterations = 50;
    options.function_tolerance = 1e-5;
    options.gradient_tolerance = 1e-6;
    options.parameter_tolerance = 1e-4;
    if (!debug) options.logging_type = ceres::SILENT;
    ceres::Solve(options, &problem, &summary);
    if (debug || !summary.IsSolutionUsable()) {
        std::cout << summary.FullReport() << std::endl;
    }
    residual = summary.final_cost;
    p[0] *= scale;
    if (debug) {
        std::cout << "Finished Direct Algo, p: " << p.transpose() << std::endl;
    }
    return p;
};

} // End namespace QI