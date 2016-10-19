/*
 *  CEST.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <unsupported/Eigen/Splines>
#include "QI/Algorithms/CEST.h"
#include "cppoptlib/solver/lbfgsbsolver.h"

namespace QI {

Eigen::ArrayXd Lorentzian(const double A, const double f0, const double w, const Eigen::ArrayXd &f) {
    Eigen::ArrayXd s = 1. - A/(1. + (2.*(f - f0)/w).square());
    return s;
}

class ZCost : public cppoptlib::BoundedProblem<double, 3> {
public:
    using BoundedProblem<double, 3>::BoundedProblem;
    using typename Problem<double, 3>::TVector;
    Eigen::ArrayXd m_frqs, m_zspec;

    ZCost(const Eigen::ArrayXd &f, const Eigen::ArrayXd &z,
          const Eigen::ArrayXd &l, const Eigen::ArrayXd &u) :
          BoundedProblem(l, u),
          m_frqs(f), m_zspec(z)
    {}

    Eigen::ArrayXd residuals(const Eigen::VectorXd &p) const {
        Eigen::ArrayXd s = Lorentzian(p[0], p[1], p[2], m_frqs);
        Eigen::ArrayXd diff = (s - m_zspec);
        Eigen::IOFormat sFormat(3);
        return diff;
    }

    double value(const TVector &p) {
        return residuals(p).square().sum();
    }
};


void CESTAlgo::apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                     std::vector<TOutput> &outputs, TConst &residual,
                     TInput &resids, TIters &its) const
{
    size_t full = m_half*2+1;
    const Eigen::Map<const Eigen::ArrayXf> z_spec(inputs[0].GetDataPointer(), full);

    auto stop = cppoptlib::Criteria<double>::defaults();
    stop.iterations = 100;
    stop.gradNorm = 1e-8;
    cppoptlib::LbfgsbSolver<ZCost> solver;
    solver.setStopCriteria(stop);
    Eigen::Array4d lower; lower << 0.0, m_ifrqs.minCoeff(), 0.001;
    Eigen::Array4d upper; upper << 1.0, m_ifrqs.maxCoeff(), (m_ifrqs.maxCoeff() - m_ifrqs.minCoeff());
    ZCost cost(m_ifrqs.segment(10,20).cast<double>(), z_spec.segment(10,20).cast<double>() / z_spec.maxCoeff(), lower, upper);

    Eigen::Vector3d p; p << 0.5, 0.0, 1.0;
    //std::cout << "p start: " << p.transpose() << std::endl;
    solver.minimize(cost, p);
    //std::cout << "p end: " << p.transpose() << std::endl;
    //std::cout << "Status " << solver.status() << " current " << solver.criteria().iterations << std::endl << " stop " << stop << std::endl;
    const float f0 = p[1];
    //std::cout << "Min f0: " << f0 << " fitted f0: " << p[1] << std::endl;
    outputs.at(2)[0] = f0;

    typedef Eigen::Spline<float, 1> TSpline;
    typedef Eigen::SplineFitting<TSpline> TFit;
    const float maxfrq = m_ifrqs.maxCoeff();
    const float minfrq = m_ifrqs.minCoeff();
    const float w = maxfrq - minfrq;
    const Eigen::ArrayXf scaledfrqs = (m_ifrqs - minfrq) / w;
    Eigen::DenseIndex degree = std::min<int>(m_ifrqs.rows() - 1, 3);
    TSpline spline = TFit::Interpolate(z_spec.transpose(), degree, scaledfrqs);
    const float ref = z_spec.maxCoeff();
    for (int f = 0; f < m_ofrqs.rows(); f++) {
        const float frq = (m_ofrqs.coeffRef(f) + f0 - minfrq)/w;
        const TSpline::PointType val = spline(frq);
        outputs.at(0)[f] = val[0];
    }
    for (int f = 0; f < m_afrqs.rows(); f++) {
        const float pfrq = (f0 + m_afrqs.coeffRef(f) - minfrq)/w;
        const float nfrq = (f0 - m_afrqs.coeffRef(f) - minfrq)/w;
        const float pos = spline(pfrq)[0];
        const float neg = spline(nfrq)[0];
        outputs.at(1)[f] = ((pos - neg)/ref);
    }
}

} // End namespace QI