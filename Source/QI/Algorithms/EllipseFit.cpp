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
#include "cppoptlib/solver/lbfgsbsolver.h"
#include "cppoptlib/solver/neldermeadsolver.h"

namespace QI {
class EllipseCost : public cppoptlib::BoundedProblem<double, 4> {
public:
    using BoundedProblem<double, 4>::BoundedProblem;
    using typename Problem<double, 4>::TVector;
    
    Eigen::ArrayXcd m_data;
    Eigen::ArrayXd m_flip, m_pincs;
    double m_TR;

    Eigen::ArrayXd residuals(const Eigen::VectorXd &p) const {
        Eigen::ArrayXcd s = QI::One_SSFP_Echo(m_flip, m_pincs, m_TR, p[0], p[1], p[2], p[3]/m_TR, 1.0);
        Eigen::ArrayXd diff = (s - m_data).abs();
        Eigen::IOFormat sFormat(3);
        std::cout << "o " << m_data.transpose().format(sFormat) << std::endl;
        std::cout << "s " << s.transpose().format(sFormat) << std::endl;
        std::cout << "d " << diff.transpose().format(sFormat) << std::endl;
        return diff;
    }

    double value(const TVector &p) {
        return residuals(p).square().sum();
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
    auto stop = cppoptlib::Criteria<double>::defaults();
    stop.iterations = 100;
    cppoptlib::NelderMeadSolver<EllipseCost> solver;
    solver.setStopCriteria(stop);
    EllipseCost cost;
    cost.m_data = data;
    cost.m_pincs = m_sequence->phase_incs();
    cost.m_flip = Eigen::ArrayXd::Constant(m_sequence->phase_incs().rows(), flip);
    cost.m_TR = m_sequence->TR();
    Eigen::Array4d lower; lower << 1.*abs(gs), 0.1, 0.001, 0.9 * arg(gs) / M_PI;
    Eigen::Array4d upper; upper << 20.*abs(gs), 4.3, 2.000, 1.1 * arg(gs) / M_PI;
    cost.setLowerBound(lower);
    cost.setUpperBound(upper);
    Eigen::Vector4d p; p << 10.*abs(gs), 1.0, 0.1, arg(gs) / M_PI;
    std::cout << "p start: " << p.transpose() << std::endl;
    solver.minimize(cost, p);
    std::cout << "p end: " << p.transpose() << std::endl;
    std::cout << "Status " << solver.status() << " current " << solver.criteria().iterations << std::endl << " stop " << stop << std::endl;
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