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

#include <array>
#include <complex>
#include <memory>
#include <string>
#include <vector>

#include "Eigen/Core"

#include "FitFunction.h"
#include "Macro.h"
#include "Model.h"
#include "SSFPSequence.h"
#include "Util.h"

using namespace std::literals;

namespace QI {

struct EllipseModel {
    using SequenceType  = QI::SSFPSequence;
    using DataType      = std::complex<double>;
    using ParameterType = double;

    static constexpr int NV = 5;
    static constexpr int ND = 0;
    static constexpr int NF = 0;

    using VaryingArray = QI_ARRAYN(ParameterType, NV);
    using FixedArray   = QI_ARRAYN(ParameterType, NF);
    static std::array<const std::string, NV> varying_names;
    static std::array<const std::string, NF> fixed_names;
    static const FixedArray                  fixed_defaults;

    const SequenceType &sequence;

    EllipseModel(const SequenceType &s) : sequence{s} {}

    // QI_ARRAYN(double, NV) bounds_lo = QI_ARRAYN(double,
    // NV)::Constant(-std::numeric_limits<double>::infinity()); QI_ARRAYN(double, NV) bounds_hi =
    // QI_ARRAYN(double, NV)::Constant(std::numeric_limits<double>::infinity());

    auto signal(const VaryingArray &v, const FixedArray & /*Unused*/) const
        -> QI_ARRAY(std::complex<double>);

    template <typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &v, const FixedArray & /*Unused*/) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T               = typename Derived::Scalar;
        using ArrayXT         = Eigen::Array<T, Eigen::Dynamic, 1>;
        const T &     G       = v[0];
        const T &     a       = v[1];
        const T &     b       = v[2];
        const T &     theta0  = v[3];
        const T &     psi0    = v[4];
        const ArrayXT theta   = theta0 - sequence.PhaseInc;
        const T       psi     = theta0 / 2.0 + psi0;
        const ArrayXT cos_th  = cos(theta);
        const ArrayXT sin_th  = sin(theta);
        const T       cos_psi = cos(psi);
        const T       sin_psi = sin(psi);
        const ArrayXT re_m =
            (cos_psi - a * (cos_th * cos_psi - sin_th * sin_psi)) * G / (1.0 - b * cos_th);
        const ArrayXT im_m =
            (sin_psi - a * (cos_th * sin_psi + sin_th * cos_psi)) * G / (1.0 - b * cos_th);
        ArrayXT result(re_m.rows() + im_m.rows());
        result << re_m, im_m;
        return result;
    }
};

using EllipseFit = BlockFitFunction<EllipseModel>;

} // End namespace QI

#endif // QI_ELLIPSE_ALGO_H