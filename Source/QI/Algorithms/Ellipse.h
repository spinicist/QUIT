/*
 *  Ellipse.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_ELLIPSE_H
#define QI_ELLIPSE_H

#include <memory>
#include <complex>
#include <array>
#include <vector>
#include <string>

#include "Eigen/Dense"

#include "QI/Macro.h"
#include "QI/Types.h"
#include "QI/Util.h"
#include "QI/Sequences/SteadyStateSequence.h"

namespace QI {

class ESAlgo : public QI::ApplyVectorXFVectorF::Algorithm {
protected:
    bool m_phaseFirst = false, m_debug = false;
    std::shared_ptr<QI::SSFPEcho> m_sequence = nullptr;
    TOutput m_zero;
public:
    typedef Eigen::Matrix<double, 6, 6> Matrix6d;
    typedef Eigen::Matrix<double, 6, 1> Vector6d;

    ESAlgo(std::shared_ptr<QI::SSFPEcho> &seq, bool debug, bool phase) :
        m_sequence(seq), m_debug(debug), m_phaseFirst(phase)
    {
        m_zero = TOutput(m_sequence->flip().rows());
        m_zero.Fill(0.);
    }

    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 1; }
    size_t numOutputs() const override { return 6; }
    size_t dataSize() const override { return m_sequence->size(); }
    size_t outputSize(const int i) const override { return m_sequence->flip().rows(); }
    void setReorderPhase(const bool p) { m_phaseFirst = p; }
    virtual std::vector<float> defaultConsts() const override {
        std::vector<float> def(1, 1.0f); // B1
        return def;
    }
    virtual const TOutput &zero(const size_t i) const override { return m_zero; }
    const std::vector<std::string> & names() const {
        static std::vector<std::string> _names = {"M", "T1", "T2", "f0", "a", "b"};
        return _names;
    }
    virtual std::array<float, 6> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &indata, const double TR, const double flip) const = 0;
    virtual bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                       std::vector<TOutput> &outputs, TConst &residual,
                       TInput &resids, TIters &its) const override
    {
        const double B1 = consts[0];
        size_t phase_stride = m_sequence->flip().rows();
        size_t flip_stride = 1;
        if (m_phaseFirst)
            std::swap(phase_stride, flip_stride);
        for (int f = 0; f < m_sequence->flip().rows(); f++) {
            const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> vf(inputs[0].GetDataPointer() + f*flip_stride, m_sequence->phase_incs().rows(), Eigen::InnerStride<>(phase_stride));
            if (m_debug) {
                std::cout << "Flip: " << m_sequence->flip() << " B1: " << B1 << " B1*flip: " << B1*m_sequence->flip() << std::endl;
            }
            std::array<float, 6> tempOutputs = this->applyFlip(vf, m_sequence->TR(), B1 * m_sequence->flip()[f]);
            for (int o = 0; o < 6; o++) {
                outputs[o][f] = tempOutputs[o];
            }
        }
        return true;
    }
};

class HyperEllipse : public ESAlgo {
public:
    HyperEllipse(std::shared_ptr<QI::SSFPEcho> &seq, bool debug, bool phase) :
        ESAlgo(seq, debug, phase)
    {}

protected:
    Eigen::MatrixXd buildS(const Eigen::ArrayXd &x, const Eigen::ArrayXd &y) const;    
    Eigen::MatrixXd fitzC() const;
    Eigen::MatrixXd hyperC(const Eigen::ArrayXd &x, const Eigen::ArrayXd &y) const;

    virtual std::array<float, 6> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &indata,
                                           const double TR, const double flip) const;
};

class ConstrainedEllipse : public ESAlgo {
public:
    ConstrainedEllipse(std::shared_ptr<QI::SSFPEcho> &seq, bool debug, bool phase, bool block) :
        ESAlgo(seq, debug, phase), m_reorderBlock(block)
    {}

protected:
    bool m_reorderBlock;

    virtual void fit(const Eigen::ArrayXd &x, const Eigen::ArrayXd &y, const Eigen::Vector2d p, const Eigen::Vector2d q, 
                     Eigen::Vector2d gammaBound, Eigen::Matrix2d &A, Eigen::Vector2d &x_c, double &g) const;

    virtual std::array<float, 6> applyFlip(const Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>> &indata,
                                           const double TR, const double flip) const;
};

} // End namespace QI

#endif // QI_ELLIPSE_H