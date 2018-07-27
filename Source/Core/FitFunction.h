/*
 *  FitFunction.h - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_FITFUNCTION_H
#define QI_FITFUNCTION_H

#include <tuple>
#include <string>
#include <Eigen/Core>
#include <itkIndex.h>
#include "ceres/ceres.h"
#include "Macro.h"
#include "Model.h"

namespace QI {

/*
 *  Return type required by the fit() function objects
 */
using FitReturnType = std::tuple<bool, std::string>;

template <typename Model_,
          bool Blocked_ = false,
          bool Indexed_ = false>
struct FitFunctionBase {
    using ModelType    = Model_;
    using SequenceType = typename ModelType::SequenceType;
    static const bool Blocked = Blocked_;
    static const bool Indexed = Indexed_;
    const SequenceType *sequence;
    ModelType     model;
    int max_iterations = 30;
    int n_inputs() const  { return 1; }
    int input_size(const int /* Unused */) const { return sequence->size(); }
    int n_fixed() const { return ModelType::NF; }
    int n_outputs() const { return ModelType::NV; }

    FitFunctionBase() = default;
    FitFunctionBase(SequenceType *s) { sequence = s; }
};

template <typename ModelType, typename FlagType_ = int>
struct FitFunction : FitFunctionBase<ModelType, false, false> {
    using Super = FitFunctionBase<ModelType, false, false>;
    using Super::Super;
    using InputType    = typename ModelType::DataType;
    using OutputType   = typename ModelType::ParameterType ;
    using ResidualType = typename ModelType::DataType;
    using FlagType     = FlagType_;   // Iterations

    virtual FitReturnType fit(const std::vector<QI_ARRAY(InputType)> &inputs,
                              const Eigen::ArrayXd &fixed, QI_ARRAYN(OutputType, ModelType::NV) &outputs,
                              ResidualType &residual, std::vector<QI_ARRAY(ResidualType)> &point_residuals, FlagType &flag) const = 0;
};

template <typename ModelType_, typename FT = int>
struct ScaledNLLSFitFunction : FitFunction<ModelType_, FT> {
    int max_iterations = 30;
    FitReturnType fit(const std::vector<QI_ARRAY(typename ScaledNLLSFitFunction::InputType)> &inputs,
                      const Eigen::ArrayXd &fixed, QI_ARRAYN(typename ScaledNLLSFitFunction::OutputType, ScaledNLLSFitFunction::ModelType::NV) &p,
                      typename ScaledNLLSFitFunction::ResidualType &residual,
                      std::vector<QI_ARRAY(typename ScaledNLLSFitFunction::ResidualType)> &residuals,
                      typename ScaledNLLSFitFunction::FlagType &iterations) const override
    {
        const double &scale = inputs[0].maxCoeff();
        if (scale < std::numeric_limits<double>::epsilon()) {
            p = ScaledNLLSFitFunction::ModelType::VaryingArray::Zero();
            residual = 0;
            return std::make_tuple(false, "Maximum data value was not positive");
        }
        const Eigen::ArrayXd data = inputs[0] / scale;
        p << 10., 1., 0.1, 12.e-6, 5.0, 0.1;
        ceres::Problem problem;
        using Cost     = ModelCost<typename ScaledNLLSFitFunction::ModelType>;
        using AutoCost = ceres::AutoDiffCostFunction<Cost, ceres::DYNAMIC, ScaledNLLSFitFunction::ModelType::NV>;
        auto *cost = new Cost(this->model, this->sequence, fixed, data);
        auto *auto_cost = new AutoCost(cost, this->sequence->size());
        problem.AddResidualBlock(auto_cost, NULL, p.data());
        problem.SetParameterLowerBound(p.data(), 0, this->model.bounds_lo[0] / scale);
        problem.SetParameterUpperBound(p.data(), 0, this->model.bounds_hi[0] / scale);
        for (int i = 1; i < ScaledNLLSFitFunction::ModelType::NV; i++) {
            problem.SetParameterLowerBound(p.data(), i, this->model.bounds_lo[i]);
            problem.SetParameterUpperBound(p.data(), i, this->model.bounds_hi[i]);
        }
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations = this->max_iterations;
        options.function_tolerance = 1e-5;
        options.gradient_tolerance = 1e-6;
        options.parameter_tolerance = 1e-4;
        options.logging_type = ceres::SILENT;
        ceres::Solve(options, &problem, &summary);
        p[0] = p[0] * scale;
        if (!summary.IsSolutionUsable()) {
            return std::make_tuple(false, summary.FullReport());
        }
        iterations = summary.iterations.size();
        residual = summary.final_cost * scale;
        if (residuals.size() > 0) {
            std::vector<double> r_temp(data.size());
            problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
            for (size_t i = 0; i < r_temp.size(); i++)
                residuals[0][i] = r_temp[i] * scale;
        }
        return std::make_tuple(true, "");
    }
};

template <typename ModelType, typename FlagType_ = int>
struct BlockFitFunction : FitFunctionBase<ModelType, true, false> {
    using Super = FitFunctionBase<ModelType, true, false>;
    using Super::Super;
    using InputType    = typename ModelType::DataType;
    using OutputType   = typename ModelType::ParameterType ;
    using ResidualType = typename ModelType::DataType;
    using FlagType     = FlagType_;   // Iterations
    using SequenceType = typename ModelType::SequenceType;
    virtual FitReturnType fit(const std::vector<QI_ARRAY(InputType)> &inputs,
                              const Eigen::ArrayXd &fixed, QI_ARRAYN(OutputType, ModelType::NV) &outputs,
                              ResidualType &residual, std::vector<QI_ARRAY(ResidualType)> &point_residuals, FlagType &flag,
                              const int block) const = 0;
};

template <typename ModelType, typename FlagType_ = int>
struct IndexedFitFunction : FitFunctionBase<ModelType, false, true> {
    using Super = FitFunctionBase<ModelType, false, true>;
    using Super::Super;
    using InputType    = typename ModelType::DataType;
    using OutputType   = typename ModelType::ParameterType ;
    using ResidualType = typename ModelType::DataType;
    using FlagType     = FlagType_;   // Iterations
    using SequenceType = typename ModelType::SequenceType;
    virtual FitReturnType fit(const std::vector<QI_ARRAY(InputType)> &inputs,
                              const Eigen::ArrayXd &fixed, QI_ARRAYN(OutputType, ModelType::NV) &outputs,
                              ResidualType &residual, std::vector<QI_ARRAY(ResidualType)> &point_residuals, FlagType &flag,
                              const itk::Index<3> &index) const = 0;
};

} // End namespace QI

#endif // QI_FITFUNCTION_H