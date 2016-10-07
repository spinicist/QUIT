/*
 *  CEST.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_CEST_H
#define QI_CEST_H

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

class CESTAlgo : public QI::ApplyVectorF::Algorithm {
protected:
    Eigen::ArrayXf m_ifrqs, m_ofrqs;
    size_t m_half;
    TOutput m_zero;
public:
    CESTAlgo(const Eigen::ArrayXf &ifrqs, const Eigen::ArrayXf &ofrqs) :
        m_ifrqs(ifrqs),
        m_ofrqs(ofrqs)
    {
        m_half = (m_ifrqs.rows() - 1) / 2;
        m_zero = TOutput(m_half);
        m_zero.Fill(0.);
    }
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 0; }
    size_t numOutputs() const override { return 1; }
    size_t dataSize() const override { return m_ifrqs.rows(); }
    size_t outputSize(const int i) const override { return m_half; }
    virtual std::vector<float> defaultConsts() const override {
        std::vector<float> def(0, 1.0f);
        return def;
    }
    virtual const TOutput &zero(const size_t i) const override { return m_zero; }
    const std::vector<std::string> & names() const {
        static std::vector<std::string> _names = {"asym"};
        return _names;
    }
    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                       std::vector<TOutput> &outputs, TConst &residual,
                       TInput &resids, TIters &its) const override;
};

} // End namespace QI

#endif // QI_CEST_H