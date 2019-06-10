/*
 *  Model.h - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_MODEL_H
#define QI_MODEL_H

#include "ImageTypes.h"
#include "Macro.h"
#include <array>
#include <string>

namespace QI {

/*
 *  Standard Model interface
 */
template <int NV_, int NF_, typename Sequence_, typename DataType_ = double> struct Model {
    using SequenceType  = Sequence_;
    using DataType      = DataType_;
    using ParameterType = double;

    static constexpr int NV = NV_; // Number of varying parameters
    static constexpr int NF = NF_; // Number of fixed parameters (fixed per voxel, e.g. B1)
    static constexpr int ND = 0;   // Number of derived parameters (calculated from varying)
    static constexpr int NI = 1;   // Number of inputs
    using VaryingArray      = QI_ARRAYN(ParameterType, NV);
    using FixedArray        = QI_ARRAYN(ParameterType, NF);
    static std::array<const std::string, NV> varying_names;
    static std::array<const std::string, NF> fixed_names;
    static const FixedArray                  fixed_defaults;

    VaryingArray bounds_lo = VaryingArray::Constant(
        1.e-12); // Don't actually use zero because that causes problems for Ceres
    VaryingArray bounds_hi = VaryingArray::Constant(std::numeric_limits<ParameterType>::infinity());

    const SequenceType &sequence;

    Model(const SequenceType &s) : sequence{s} {}

    template <typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &varying, const FixedArray &fixed) const
        -> QI_ARRAY(typename Derived::Scalar);
};

/*
 *  A generic Ceres Cost Function compatible with auto-differentation
 */
template <typename Model> struct ModelCost {
    using Sequence     = typename Model::SequenceType;
    using VaryingArray = typename Model::VaryingArray;
    using FixedArray   = typename Model::FixedArray;
    using DataArray    = QI_ARRAY(typename Model::DataType);
    const Model &    model;
    const FixedArray fixed;
    const QI_ARRAY(typename Model::DataType) data;

    ModelCost(const Model &m, const FixedArray &f, const DataArray &d) :
        model(m), fixed(f), data(d) {}

    template <typename T> bool operator()(const T *const vin, T *rin) const {
        Eigen::Map<QI_ARRAY(T)>                         r(rin, data.rows());
        const Eigen::Map<const QI_ARRAYN(T, Model::NV)> v(vin);
        const auto                                      calc = model.signal(v, fixed);
        r                                                    = data - calc;
        return true;
    }
};

} // End namespace QI

#endif // MODEL_H