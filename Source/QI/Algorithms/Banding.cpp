/*
 *  Banding.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Algorithms/Banding.h"

namespace QI {

void BandAlgo::setPhases(const size_t p) {
    if (p < 4)
        QI_EXCEPTION("Must have a minimum of 4 phase-cycling patterns.");
    if ((p % 2) != 0)
        QI_EXCEPTION("Number of phases must be even.");
    m_phases = p;
    m_lines = m_phases / 2;;
}

void BandAlgo::setInputSize(const size_t s) {
    m_flips = s / m_phases;
    m_zero.SetSize(m_flips);
    m_zero.Fill(std::complex<float>(0.));
}

std::vector<float> BandAlgo::defaultConsts() const {
    std::vector<float> def;
    return def;
}
const BandAlgo::TOutput &BandAlgo::zero(const size_t i) const {
    return m_zero;
}

void BandAlgo::apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                     std::vector<TOutput> &outputs, TConst &residual,
                     TInput &resids, TIters &its) const
{
    size_t phase_stride = m_flips;
    size_t flip_stride = 1;
    if (m_phaseFirst)
        std::swap(phase_stride, flip_stride);
    for (int f = 0; f < m_flips; f++) {
        const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> vf(inputs[0].GetDataPointer() + f*flip_stride, m_phases, Eigen::InnerStride<>(phase_stride));
        outputs[0][f] = this->applyFlip(vf);
    }
}

std::complex<float> GeometricSolution(const Eigen::ArrayXcd &a, const Eigen::ArrayXcd &b, RegEnum regularise) {
    eigen_assert(a.rows() == b.rows());
    std::complex<double> sum(0., 0.);
    double N = 0;
    for (size_t i = 0; i < a.rows(); i++) {
        for (size_t j = i + 1; j < a.rows(); j++) {
            const std::complex<double> di = b[i] -  a[i], dj = b[j] - a[j];
            const std::complex<double> ni(-di.imag(), di.real()), nj(-dj.imag(), dj.real());

            const double mu = QI::cdot(a[j] - a[i], nj) / QI::cdot(di, nj);
            const double nu = QI::cdot(a[i] - a[j], ni) / QI::cdot(dj, ni);
            const double xi = 1.0 - pow(QI::cdot(di, dj) / (abs(di)*abs(dj)), 2.0);

            const std::complex<double> cs = (a[i] + a[j] + b[i] + b[j]) / 4.0;
            const std::complex<double> gs = a[i] + mu * di;

            switch (regularise) {
            case RegEnum::None: sum += gs; break;
            case RegEnum::Magnitude:
                if (norm(gs) < std::max({std::norm(a[i]), std::norm(a[j]), std::norm(b[i]), std::norm(b[j])})) {
                    sum += gs;
                } else {
                    sum += cs;
                }
                break;
            case RegEnum::Line:
                if ((mu > -xi) && (mu < 1 + xi) && (nu > -xi) && (nu < 1 + xi)) {
                    sum += gs;
                } else {
                    sum += cs;
                }
                break;
            }
            N += 1;
        }
    }
    return static_cast<std::complex<float>>(sum / N);
}

std::complex<float> GSAlgo::applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const {
    Eigen::ArrayXcd a(m_lines);
    Eigen::ArrayXcd b(m_lines);
    Eigen::ArrayXcd full = vf.cast<std::complex<double>>();
    SplitBlocks(full, a, b, m_reorderBlock);
    return GeometricSolution(a, b, m_Regularise);
}

} // End namespace QI