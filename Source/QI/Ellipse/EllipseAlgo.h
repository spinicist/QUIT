/*
 *  EllipseAlgo.h
 *
 *  Copyright (c) 2016, 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_ELLIPSE_ALGO_H
#define QI_ELLIPSE_ALGO_H

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
#include "Filters/ApplyAlgorithmFilter.h"

namespace QI {

class EllipseAlgo : public QI::ApplyVectorXFVectorF::Algorithm {
protected:
    bool m_debug = false;
    std::shared_ptr<QI::SSFPEcho> m_sequence = nullptr;
    TOutput m_zero;
    virtual Eigen::ArrayXd apply_internal(const Eigen::ArrayXcf &input, const double flip, const double TR, const Eigen::ArrayXd &phi, const bool debug, float &residual) const = 0;
public:
    EllipseAlgo(std::shared_ptr<QI::SSFPEcho> &seq, bool debug);

    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 1; }
    size_t dataSize() const override { return m_sequence->size(); }
    size_t outputSize(const int i) const override { return m_sequence->flip().rows(); }
    virtual std::vector<float> defaultConsts() const override {
        std::vector<float> def(1, 1.);
        return def;
    }
    virtual const TOutput &zero(const size_t i) const override { return m_zero; }
    virtual bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                       std::vector<TOutput> &outputs, TConst &residual,
                       TInput &resids, TIters &its) const override;
    
    virtual const std::vector<std::string> & names() const = 0;
};

} // End namespace QI

#endif // QI_ELLIPSE_ALGO_H