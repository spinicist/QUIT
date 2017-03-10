/*
 *  Ellipse.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Algorithms/EllipseFit.h"
#include "QI/Algorithms/Banding.h"
#include "ceres/ceres.h"

namespace QI {
class EllipseCost {
private:
    const Eigen::ArrayXcd &m_data;
    std::shared_ptr<QI::SSFPEcho> m_s;
public:
    EllipseCost(const Eigen::ArrayXcd &d, std::shared_ptr<QI::SSFPEcho> s) :
        m_data(d), m_s(s)
    {
    }

    bool operator() (double const* const* p, double* resids) const {
        Eigen::ArrayXcd s = QI::One_SSFP_Echo(m_s->allFlip(), m_s->allPhi(), m_s->TR(), p[0][0], p[0][1], p[0][2], p[0][3]/m_s->TR(), 1.0);
        Eigen::Map<Eigen::ArrayXd> r(resids, m_data.size());
        r = (s - m_data).abs();
        return true;
    }
};

std::array<float, 6> FitEllipse::applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &indata,
                                           const double TR, const double flip) const {
    Eigen::ArrayXcd data = indata.cast<std::complex<double>>();
    const double scale = data.abs().maxCoeff();
    data /= scale;

    Eigen::ArrayXcd a(data.rows() / 2);
    Eigen::ArrayXcd b(data.rows() / 2);
    QI::SplitBlocks(data, a, b, m_reorderBlock);
    const std::complex<double> gs = QI::GeometricSolution(a, b, QI::RegEnum::Line);

    auto *cost = new ceres::DynamicNumericDiffCostFunction<EllipseCost>(new EllipseCost(data, m_sequence));
    cost->AddParameterBlock(4);
    cost->SetNumResiduals(data.size());
    Eigen::Array4d p {10.*abs(gs), 1.0, 0.1, arg(gs) / M_PI};
    ceres::Problem problem;
    problem.AddResidualBlock(cost, NULL, p.data());
    problem.SetParameterLowerBound(p.data(), 0, 1.*abs(gs));
    problem.SetParameterLowerBound(p.data(), 0, 20.*abs(gs));
    problem.SetParameterLowerBound(p.data(), 1, 0.1);
    problem.SetParameterUpperBound(p.data(), 1, 4.3);
    problem.SetParameterLowerBound(p.data(), 2, 0.001);
    problem.SetParameterUpperBound(p.data(), 2, 4.3);
    problem.SetParameterLowerBound(p.data(), 3, 0.9 * arg(gs) / M_PI);
    problem.SetParameterUpperBound(p.data(), 3, 1.1 * arg(gs) / M_PI);
    ceres::Solver::Options options;
    ceres::Solver::Summary summary;
    options.max_num_iterations = 50;
    options.function_tolerance = 1e-5;
    options.gradient_tolerance = 1e-6;
    options.parameter_tolerance = 1e-4;
    ceres::Solve(options, &problem, &summary);
    std::array<float, 6> outputs;
    outputs[0] = p[0] * scale;
    outputs[1] = p[1];
    outputs[2] = p[2];
    outputs[3] = p[3] / m_sequence->TR();
    outputs[4] = 0;
    outputs[5] = 0;
    return outputs;
};

} // End namespace QI