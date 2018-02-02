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

#include "Macro.h"
#include "ApplyTypes.h"
#include "Util.h"
#include "SteadyStateSequence.h"

namespace QI {

class CESTAlgo : public QI::ApplyVectorF::Algorithm {
protected:
    Eigen::ArrayXf m_ifrqs, m_ofrqs, m_afrqs;
    size_t m_half;
    TOutput m_zero1, m_zero2, m_zero3;
public:
    CESTAlgo(const Eigen::ArrayXf &ifrqs,
             const Eigen::ArrayXf &ofrqs,
             const Eigen::ArrayXf &afrqs) :
        m_ifrqs(ifrqs),
        m_ofrqs(ofrqs),
        m_afrqs(afrqs)
    {
        m_half = (m_ifrqs.rows() - 1) / 2;
        m_zero1 = TOutput(m_ofrqs.rows()); m_zero2.Fill(0.);
        m_zero2 = TOutput(m_afrqs.rows()); m_zero3.Fill(0.);
        m_zero3 = TOutput(1); m_zero3.Fill(0.);
    }
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 0; }
    size_t numOutputs() const override { return 3; }
    size_t dataSize() const override { return m_ifrqs.rows(); }
    size_t outputSize(const int i) const override {
        switch (i) {
            case 0: return m_ofrqs.rows();
            case 1: return m_afrqs.rows();
            case 2: return 1;
            default: QI_EXCEPTION("Requested invalid output " << i);
        }
    }
    std::vector<float> defaultConsts() const override {
        std::vector<float> def(0, 1.0f);
        return def;
    }
    const TOutput &zero(const size_t i) const override {
        switch (i) {
            case 0: return m_zero1;
            case 1: return m_zero2;
            case 2: return m_zero3;
            default: QI_EXCEPTION("Requested invalid output " << i);
        }
    }
    const std::vector<std::string> & names() const {
        static std::vector<std::string> _names = {"centered","asym", "f0"};
        return _names;
    }
    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIterations &its) const override;
};

} // End namespace QI

#endif // QI_CEST_H