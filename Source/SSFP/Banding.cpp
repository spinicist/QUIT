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

#include "Banding.h"
#include "Log.h"

namespace QI {

/*
 * Helper Functions
 */
template <typename T> inline T cdot(const std::complex<T> &a, const std::complex<T> &b) {
    return real(a) * real(b) + imag(a) * imag(b);
}

BandFunctor::BandFunctor(const int isz, const int p, const bool rp, const bool rb)
    : m_phaseFirst(rp), m_reorderBlock(rb) {
    if (p < 4)
        QI::Fail("Must have a minimum of 4 phase-cycling patterns.");
    if ((p % 2) != 0)
        QI::Fail("Number of phases must be even.");
    m_phases = p;
    m_lines  = m_phases / 2;
    m_flips  = isz / m_phases;
}

itk::VariableLengthVector<std::complex<float>> BandFunctor::
                                               operator()(const itk::VariableLengthVector<std::complex<float>> &vec) const {
    size_t phase_stride = m_flips;
    size_t flip_stride  = 1;
    if (m_phaseFirst)
        std::swap(phase_stride, flip_stride);
    itk::VariableLengthVector<std::complex<float>> output(m_flips);
    for (size_t f = 0; f < m_flips; f++) {
        const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> vf(
            vec.GetDataPointer() + f * flip_stride, m_phases, Eigen::InnerStride<>(phase_stride));
        output[f] = this->applyFlip(vf);
    }
    return output;
}

std::complex<float> GeometricSolution(const Eigen::ArrayXcd &a, const Eigen::ArrayXcd &b,
                                      RegEnum regularise) {
    eigen_assert(a.rows() == b.rows());
    std::complex<double> sum(0., 0.);
    double               N = 0;
    for (int i = 0; i < a.rows(); i++) {
        for (int j = i + 1; j < a.rows(); j++) {
            const std::complex<double> di = b[i] - a[i], dj = b[j] - a[j];
            const std::complex<double> ni(-di.imag(), di.real()), nj(-dj.imag(), dj.real());

            const double mu = QI::cdot(a[j] - a[i], nj) / QI::cdot(di, nj);
            const double nu = QI::cdot(a[i] - a[j], ni) / QI::cdot(dj, ni);
            const double xi = 1.0 - pow(QI::cdot(di, dj) / (abs(di) * abs(dj)), 2.0);

            const std::complex<double> cs = (a[i] + a[j] + b[i] + b[j]) / 4.0;
            const std::complex<double> gs = a[i] + mu * di;

            switch (regularise) {
            case RegEnum::None:
                sum += gs;
                break;
            case RegEnum::Magnitude:
                if (norm(gs) < std::max({std::norm(a[i]), std::norm(a[j]), std::norm(b[i]),
                                         std::norm(b[j])})) {
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

std::complex<float>
CSFunctor::applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const {
    return vf.mean();
}

std::complex<float> MagMeanFunctor::applyFlip(
    const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const {
    return vf.abs().mean();
}

std::complex<float>
RMSFunctor::applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const {
    float sum = vf.abs().square().sum();
    return std::complex<float>(sqrt(sum / vf.rows()), 0.);
}

std::complex<float>
MaxFunctor::applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const {
    std::complex<float> max = std::numeric_limits<std::complex<float>>::lowest();
    for (int i = 0; i < vf.rows(); i++) {
        if (std::abs(vf[i]) > std::abs(max))
            max = vf[i];
    }
    return max;
}

std::complex<float>
GSFunctor::applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const {
    Eigen::ArrayXcd a(m_lines);
    Eigen::ArrayXcd b(m_lines);
    Eigen::ArrayXcd full = vf.cast<std::complex<double>>();
    SplitBlocks(full, a, b, m_reorderBlock);
    return GeometricSolution(a, b, m_Regularise);
}

} // End namespace QI