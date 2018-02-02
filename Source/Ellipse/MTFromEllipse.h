/*
 *  MTFromEllipseFilter.h
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_ELLIPSE_MTFROMELLIPSE_H
#define QI_ELLIPSE_MTFROMELLIPSE_H

#include "ApplyTypes.h"

namespace QI {

class MTFromEllipse : public QI::ApplyF::Algorithm {
public:
    const static size_t NumOutputs = 5;
protected:
    const Eigen::ArrayXd &flips;
    const Eigen::ArrayXd &intB1;
    const Eigen::ArrayXd &TRs;
    const Eigen::ArrayXd &TRFs;
    const double T2r;
    const bool debug;
public:
    MTFromEllipse(const Eigen::ArrayXd &f, const Eigen::ArrayXd &iB, const Eigen::ArrayXd &tr, const Eigen::ArrayXd &trf, const double T2, const bool d);
    size_t numInputs() const override { return 3; }
    size_t numConsts() const override { return 2; }
    size_t numOutputs() const override { return NumOutputs; }
    size_t dataSize() const override { return (flips.size() * 3); }
    const std::vector<std::string> &names() const {
        static std::vector<std::string> _names = {"M0", "F", "kf", "T1f", "T2f"};
        return _names;
    }
    std::vector<float> defaultConsts() const override;
    TOutput zero() const override { return 0.f; }
    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIterations &its) const override;
};

} // End namespace QI

#endif // QI_ELLIPSE_MTFROMELLIPSE_H