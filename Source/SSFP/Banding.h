/*
 *  Banding.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_BANDING_H
#define QI_BANDING_H

#include <complex>
#include "Eigen/Core"
#include "itkVariableLengthVector.h"

namespace QI {

template<typename Derived>
void SplitBlocks(const Eigen::ArrayBase<Derived> &full, Eigen::ArrayBase<Derived> &a, Eigen::ArrayBase<Derived> &b, const bool reorder) {
    if (reorder) {
        for (int i = 0; i < a.rows(); i++) {
            a[i] = static_cast<std::complex<double>>(full[i*2]);
            b[i] = static_cast<std::complex<double>>(full[i*2+1]);
        }
    } else {
        a = full.head(a.rows());
        b = full.tail(b.rows());
    }
}

enum class RegEnum { None = 0, Line, Magnitude };
std::complex<float> GeometricSolution(const Eigen::ArrayXcd &a, const Eigen::ArrayXcd &b, RegEnum r);

struct BandFunctor {
    size_t m_flips, m_lines, m_crossings, m_phases = 4;
    bool m_phaseFirst = false, m_reorderBlock = false;

    BandFunctor() = default;
    BandFunctor(const int isz, const int p, const bool rp, const bool rb);

    bool operator!=(const BandFunctor &) const { return true; }
    bool operator==(const BandFunctor &other) const { return !(*this != other); }
    itk::VariableLengthVector<std::complex<float>> operator()(const itk::VariableLengthVector<std::complex<float>> &vec) const;
    virtual std::complex<float> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const = 0;
};

struct CSFunctor : BandFunctor {
    using BandFunctor::BandFunctor;
    std::complex<float> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const override;
};

struct MagMeanFunctor : BandFunctor {
    using BandFunctor::BandFunctor;
    std::complex<float> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const override;
};

struct RMSFunctor : BandFunctor {
    using BandFunctor::BandFunctor;
    std::complex<float> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const override;
};

struct MaxFunctor : BandFunctor {
    using BandFunctor::BandFunctor;
    std::complex<float> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const override;
};

struct GSFunctor : BandFunctor {
    using BandFunctor::BandFunctor;
    RegEnum m_Regularise = RegEnum::Line;
    const RegEnum &regularise()       { return m_Regularise; }
    void setRegularise(const RegEnum &r) { m_Regularise = r;}
    std::complex<float> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const override;
};

} // End namespace QI

#endif // QI_BANDING_H