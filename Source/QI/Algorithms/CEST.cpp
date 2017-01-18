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

    // Find closest indices to -2/+2 PPM and only fit Lorentzian between them
    Eigen::ArrayXf::Index indP2, indM2;
    (m_ifrqs + 2.0).abs().minCoeff(&indM2);
    (m_ifrqs - 2.0).abs().minCoeff(&indP2);
    //std::cout << "ind " << indM2 << " " << indP2 << std::endl;
    if (indM2 > indP2)
        std::swap(indM2, indP2);
    //std::cout << "ind " << indM2 << " " << indP2 << std::endl;
    Eigen::ArrayXf::Index sz = indP2 - indM2;
    Eigen::Array3d lower; lower << 0.0, -2.0, 0.001;
    Eigen::Array3d upper; upper << 1.0,  2.0, 4.0;

    /*std::cout << "seg " << m_ifrqs.segment(indM2, sz).transpose() << std::endl;
    std::cout << "zpc " << z_spec.segment(indM2,sz).transpose() << std::endl;
    std::cout << "lower " << lower.transpose() << " upper " << upper.transpose() << std::endl;*/
    ZCost cost(m_ifrqs.segment(indM2,sz).cast<double>(), z_spec.segment(indM2,sz).cast<double>() / z_spec.segment(indM2,sz).maxCoeff(), lower, upper);

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
    /*std::cout << "sFreqs " << scaledfrqs.transpose() << std::endl;
    std::cout << "z_spec " << z_spec.transpose() << std::endl;
    std::cout << "degree " << degree << std::endl;*/
    TSpline spline;
    if (scaledfrqs[0] > 0)
        spline = TFit::Interpolate(z_spec.reverse().transpose(), degree, scaledfrqs.reverse());
    else
        spline = TFit::Interpolate(z_spec.transpose(), degree, scaledfrqs);
    const float ref = z_spec.maxCoeff();
    for (int f = 0; f < m_ofrqs.rows(); f++) {
        const float frq = (m_ofrqs.coeffRef(f) + f0 - minfrq)/w;
        const TSpline::PointType val = spline(frq);
        outputs.at(0)[f] = val[0];
        //std::cout << "f " << f << " frq " << frq << " val " << val.transpose() << std::endl;
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