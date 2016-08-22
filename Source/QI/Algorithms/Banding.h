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
#include "Eigen/Dense"

#include "QI/Macro.h"
#include "QI/Types.h"
#include "QI/Util.h"

namespace QI {

/*
 * Helper Functions
 */
template<typename T> inline T cdot(const std::complex<T> &a, const std::complex<T> &b) {
    return real(a * conj(b));
}

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

class BandAlgo : public ApplyVectorXF::Algorithm {
protected:
    size_t m_flips, m_lines, m_crossings, m_phases = 4;
    bool m_phaseFirst = false, m_reorderBlock = false;
    TOutput m_zero;
public:
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 0; }
    size_t numOutputs() const override { return 1; }
    size_t dataSize() const override { return m_flips * m_phases; }
    size_t outputSize(const int i) const override { return m_flips; }
    void setPhases(const size_t p);
    void setInputSize(const size_t s);
    void setReorderPhase(const bool p) { m_phaseFirst = p; }
    void setReorderBlock(const bool b) { m_reorderBlock = b; }
    virtual std::vector<float> defaultConsts() const override;
    virtual const TOutput &zero(const size_t i) const override;
    virtual std::complex<float> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const = 0;
    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                       std::vector<TOutput> &outputs, TConst &residual,
                       TInput &resids, TIters &its) const override;
};

class CSAlgo : public BandAlgo {
public:
    std::complex<float> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const override {
        return vf.mean();
    }
};

class MagMeanAlgo : public BandAlgo {
public:
    std::complex<float> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const override {
        return vf.abs().mean();
    }
};

class RMSAlgo : public BandAlgo {
public:
    std::complex<float> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const override {
        float sum = vf.abs().square().sum();
        return std::complex<float>(sqrt(sum / vf.rows()), 0.);
    }
};

class MaxAlgo : public BandAlgo {
public:
    std::complex<float> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const override {
        std::complex<float> max = std::numeric_limits<std::complex<float>>::lowest();
        for (size_t i = 0; i < vf.rows(); i++) {
            if (std::abs(vf[i]) > std::abs(max)) max = vf[i];
        }
        return max;
    }
};

class GSAlgo : public BandAlgo {
protected:
    RegEnum m_Regularise = RegEnum::Line;
public:
    const RegEnum &regularise()       { return m_Regularise; }
    void setRegularise(const RegEnum &r) { m_Regularise = r;}
    std::complex<float> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const override;
};

} // End namespace QI

#endif // QI_BANDING_H