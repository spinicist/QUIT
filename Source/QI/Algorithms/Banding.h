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

template<typename T> inline T cdot(const std::complex<T> &a, const std::complex<T> &b) {
    return real(a * conj(b));
}

class BandAlgo : public ApplyVectorXF::Algorithm {
protected:
    size_t m_flips, m_lines, m_crossings, m_phases = 4;
    bool m_reorderPhase = false, m_reorderBlock = false;
    TOutput m_zero;
public:
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 0; }
    size_t numOutputs() const override { return 1; }
    size_t dataSize() const override { return m_flips * m_phases; }
    size_t outputSize(const int i) const override { return m_flips; }
    void setPhases(const size_t p) {
        if (p < 4)
            QI_EXCEPTION("Must have a minimum of 4 phase-cycling patterns.");
        if ((p % 2) != 0)
            QI_EXCEPTION("Number of phases must be even.");
        m_phases = p;
        m_lines = m_phases / 2;
        m_crossings = QI::Choose(m_lines, 2);
    }
    void setInputSize(const size_t s) {
        m_flips = s / m_phases;
        m_zero.SetSize(m_flips);
        m_zero.Fill(std::complex<float>(0.));
    }
    void setReorderPhase(const bool p) { m_reorderPhase = p; }
    void setReorderBlock(const bool b) { m_reorderBlock = b; }

    virtual std::vector<float> defaultConsts() override {
        std::vector<float> def;
        return def;
    }
    virtual const TOutput &zero(const size_t i) const override {
        return m_zero;
    }
    virtual std::complex<float> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const = 0;
    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                       std::vector<TOutput> &outputs, TConst &residual,
                       TInput &resids, TIters &its) const override
    {
        size_t phase_stride = 1;
        if (m_reorderPhase)
            phase_stride = m_flips;
        Eigen::ArrayXcf out;
        for (int f = 0; f < m_flips; f++) {
            const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> vf(inputs[0].GetDataPointer() + f, inputs[0].GetSize(), Eigen::InnerStride<>(phase_stride));
            outputs[0][f] = this->applyFlip(vf);
        }
    }
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
public:
    enum class RegEnum { None = 0, Line, Magnitude };
protected:
    RegEnum m_Regularise = RegEnum::Line;
public:
    const RegEnum &regularise()       { return m_Regularise; }
    void setRegularise(const RegEnum &r) { m_Regularise = r;}
    std::complex<float> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &vf) const override {
        Eigen::ArrayXcd a(m_flips / 2);
        Eigen::ArrayXcd b(m_flips / 2);
        if (m_reorderBlock) {
            for (int i = 0; i < m_phases; i++) {
                a[i] = static_cast<std::complex<double>>(vf[i*2]);
                b[i] = static_cast<std::complex<double>>(vf[i*2+1]);
            }
        } else {
            a = vf.head(m_lines).cast<std::complex<double>>();
            b = vf.tail(m_lines).cast<std::complex<double>>();
        }

        std::complex<double> sum(0., 0.);
        for (size_t i = 0; i < m_lines; i++) {
            for (size_t j = i + 1; j < m_lines; j++) {
                const std::complex<double> di = b[i] -  a[i], dj = b[j] - a[j];
                const std::complex<double> ni(di.imag(), -di.real()), nj(dj.imag(), -dj.real());

                const double mu = QI::cdot(a[j] - a[i], nj) / QI::cdot(di, nj);
                const double nu = QI::cdot(a[i] - a[j], ni) / QI::cdot(dj, ni);
                const double xi = 1.0 - pow(QI::cdot(di, dj) / (abs(di)*abs(dj)), 2.0);

                const std::complex<double> cs = (a[i] + a[j] + b[i] + b[j]) / 4.0;
                const std::complex<double> gs = a[i] + mu * di;

                switch (m_Regularise) {
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
            }
        }
        return static_cast<std::complex<float>>(sum / static_cast<double>(m_crossings));
    }
};

} // End namespace QI

#endif // QI_BANDING_H