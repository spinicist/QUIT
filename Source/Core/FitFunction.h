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

#pragma once

#include "Macro.h"
#include "Model.h"
#include <Eigen/Core>
#include <itkIndex.h>
#include <string>
#include <tuple>

namespace QI {

/*
 *  Return type required by the fit() function objects
 */

struct FitReturnType {
    bool        success;
    std::string message;
};

template <typename Model_, bool Blocked_ = false, bool Indexed_ = false> struct FitFunctionBase {
    using ModelType           = Model_;
    using RMSErrorType        = double;
    static const bool Blocked = Blocked_;
    static const bool Indexed = Indexed_;

    ModelType model;
    FitFunctionBase(ModelType &m) : model{m} {}

    long input_size(long const &i) const { return model.input_size(i); }
};

template <typename ModelType, typename FlagType_ = int>
struct FitFunction : FitFunctionBase<ModelType, false, false> {
    using Super = FitFunctionBase<ModelType, false, false>;
    using Super::Super;
    using typename Super::RMSErrorType;
    using InputType  = typename ModelType::DataType;
    using OutputType = typename ModelType::ParameterType;
    using FlagType   = FlagType_; // Iterations

    FitFunction(ModelType &m) : Super{m} {}

    virtual FitReturnType fit(const std::vector<QI_ARRAY(InputType)> &inputs,
                              typename ModelType::FixedArray const &  fixed,
                              typename ModelType::VaryingArray &      outputs,
                              typename ModelType::CovarArray *        cov,
                              RMSErrorType &                          rmse,
                              std::vector<QI_ARRAY(InputType)> &      residuals,
                              FlagType &                              flag) const = 0;
};

template <typename ModelType, typename FlagType_ = int>
struct NLLSFitFunction : FitFunction<ModelType> {
    using Super = FitFunction<ModelType, FlagType_>;
    using Super::Super;
    using typename Super::RMSErrorType;
    using InputType  = typename ModelType::DataType;
    using OutputType = typename ModelType::ParameterType;
    using FlagType   = FlagType_; // Iterations

    NLLSFitFunction(ModelType &m) : Super{m} {}

    FitReturnType fit(const std::vector<QI_ARRAY(InputType)> &inputs,
                      typename ModelType::FixedArray const &  fixed,
                      typename ModelType::VaryingArray &      p,
                      typename ModelType::CovarArray *        cov,
                      RMSErrorType &                          rmse,
                      std::vector<QI_ARRAY(InputType)> &      residuals,
                      FlagType &                              iterations) const {
        auto const &   data = inputs[0];
        ceres::Problem problem;
        using Cost      = ModelCost<ModelType>;
        using AutoCost  = ceres::AutoDiffCostFunction<Cost, ceres::DYNAMIC, ModelType::NV>;
        auto *cost      = new Cost{this->model, fixed, data};
        auto *auto_cost = new AutoCost(cost, this->model.sequence.size());
        problem.AddResidualBlock(auto_cost, NULL, p.data());
        for (int i = 0; i < ModelType::NV; i++) {
            problem.SetParameterLowerBound(p.data(), i, this->model.bounds_lo[i]);
            problem.SetParameterUpperBound(p.data(), i, this->model.bounds_hi[i]);
        }
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = 15;
        options.function_tolerance  = 1e-6;
        options.gradient_tolerance  = 1e-7;
        options.parameter_tolerance = 1e-5;
        options.logging_type        = ceres::SILENT;
        p << this->model.start;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            return {false, summary.FullReport()};
        }
        iterations = summary.iterations.size();

        Eigen::ArrayXd const rs  = (data - this->model.signal(p, fixed));
        double const         var = rs.square().sum();
        rmse                     = sqrt(var / data.rows());
        if (residuals.size() > 0) {
            residuals[0] = rs;
        }
        if (cov) {
            QI::GetModelCovariance<ModelType>(problem, p, var / (data.rows() - ModelType::NV), cov);
        }

        return {true, ""};
    }
};

template <typename ModelType, typename FlagType_ = int>
struct ScaledAutoDiffFit : FitFunction<ModelType, FlagType_> {
    using Super = FitFunction<ModelType, FlagType_>;
    using Super::Super;
    using typename Super::RMSErrorType;
    using InputType  = typename ModelType::DataType;
    using OutputType = typename ModelType::ParameterType;
    using FlagType   = FlagType_; // Iterations

    ScaledAutoDiffFit(ModelType &m) : Super{m} {}

    FitReturnType fit(std::vector<QI_ARRAY(InputType)> const &inputs,
                      typename ModelType::FixedArray const &  fixed,
                      typename ModelType::VaryingArray &      p,
                      typename ModelType::CovarArray *        cov,
                      RMSErrorType &                          rmse,
                      std::vector<QI_ARRAY(InputType)> &      residuals,
                      FlagType &                              iterations) const override {
        const double &scale = inputs[0].maxCoeff();
        if (scale < std::numeric_limits<double>::epsilon()) {
            p    = ModelType::VaryingArray::Zero();
            rmse = 0;
            return {false, "Maximum data value was not positive"};
        }
        const Eigen::ArrayXd data = inputs[0] / scale;
        ceres::Problem       problem;
        using Cost      = ModelCost<ModelType>;
        using AutoCost  = ceres::AutoDiffCostFunction<Cost, ceres::DYNAMIC, ModelType::NV>;
        auto *cost      = new Cost{this->model, fixed, data};
        auto *auto_cost = new AutoCost(cost, this->model.sequence.size());
        problem.AddResidualBlock(auto_cost, NULL, p.data());
        for (int i = 0; i < ModelType::NV; i++) {
            problem.SetParameterLowerBound(p.data(), i, this->model.bounds_lo[i]);
            problem.SetParameterUpperBound(p.data(), i, this->model.bounds_hi[i]);
        }
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = 15;
        options.function_tolerance  = 1e-6;
        options.gradient_tolerance  = 1e-7;
        options.parameter_tolerance = 1e-5;
        options.logging_type        = ceres::SILENT;
        p << this->model.start;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            return {false, summary.FullReport()};
        }
        iterations               = summary.iterations.size();
        Eigen::ArrayXd const rs  = (data - this->model.signal(p, fixed));
        double const         var = rs.square().sum();
        rmse                     = sqrt(var / data.rows()) * scale;
        if (residuals.size() > 0) {
            residuals[0] = rs * scale;
        }
        if (cov) {
            QI::GetModelCovariance<ModelType>(problem, p, var / (data.rows() - ModelType::NV), cov);
        }
        p[0] = p[0] * scale;
        return {true, ""};
    }
};

template <typename ModelType, int NScale = 1>
struct ScaledNumericDiffFit : FitFunction<ModelType, int> {
    using Super = FitFunction<ModelType, int>;
    using Super::Super;
    using typename Super::RMSErrorType;
    using InputType  = typename ModelType::DataType;
    using OutputType = typename ModelType::ParameterType;

    ScaledNumericDiffFit(ModelType &m) : Super{m} {}

    // This has to match the function signature that will be called in ModelFitFilter (which depends
    // on Blocked/Indexed. The return type is a simple struct indicating success, and on failure
    // also the reason for failure
    QI::FitReturnType
    fit(std::vector<Eigen::ArrayXd> const &   inputs,  // Input: signal data
        typename ModelType::FixedArray const &fixed,   // Input: Fixed parameters
        typename ModelType::VaryingArray &    varying, // Output: Varying parameters
        typename ModelType::CovarArray *      cov,
        RMSErrorType &                        rmse,      // Output: root-mean-square error
        std::vector<Eigen::ArrayXd> &         residuals, // Optional output: point residuals
        int &                                 iterations /* Usually iterations */) const override {
        // First scale down the raw data so that PD will be roughly the same magnitude as other
        // parameters This is important for numerical stability in the optimiser
        double scale = inputs[0].maxCoeff();
        if (scale < std::numeric_limits<double>::epsilon()) {
            varying = ModelType::VaryingArray::Zero();
            rmse    = 0.0;
            return {false, "Maximum data value was zero or less"};
        }
        Eigen::ArrayXd const data = inputs[0] / scale;

        // Setup Ceres
        ceres::Problem problem;
        using Diff = ceres::NumericDiffCostFunction<QI::ModelCost<ModelType>,
                                                    ceres::CENTRAL,
                                                    ceres::DYNAMIC,
                                                    ModelType::NV>;
        auto *cost = new Diff(new QI::ModelCost<ModelType>{this->model, fixed, data},
                              ceres::TAKE_OWNERSHIP,
                              this->model.sequence.size());
        auto *loss = new ceres::HuberLoss(1.0); // Don't know if this helps

        // This is where the parameters and cost functions actually get added to Ceres
        problem.AddResidualBlock(cost, loss, varying.data());

        // Set up parameter bounds
        for (int i = 0; i < ModelType::NV; i++) {
            problem.SetParameterLowerBound(varying.data(), i, this->model.lo[i]);
            problem.SetParameterUpperBound(varying.data(), i, this->model.hi[i]);
        }

        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = 15;
        options.function_tolerance  = 1e-5;
        options.gradient_tolerance  = 1e-6;
        options.parameter_tolerance = 1e-4;
        options.logging_type        = ceres::SILENT;

        varying = this->model.start;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            return {false, summary.FullReport()};
        }
        iterations = summary.iterations.size();
        double              var;
        std::vector<double> rs(data.size());
        problem.Evaluate(ceres::Problem::EvaluateOptions(), &var, &rs, nullptr, nullptr);
        rmse = sqrt(var / data.rows()) * scale;
        if (residuals.size() > 0) {
            for (long ii = 0; ii < residuals[0].size(); ii++) {
                residuals[0][ii] = rs[ii] * scale;
            }
        }
        if (cov) {
            QI::GetModelCovariance<ModelType>(
                problem, varying, var / (data.rows() - ModelType::NV), cov);
        }
        varying.template head<NScale>() *= scale; // Multiply signals/proton density back up

        return {true, ""};
    }
};

template <typename ModelType, typename FlagType_ = int>
struct BlockFitFunction : FitFunctionBase<ModelType, true, false> {
    using Super = FitFunctionBase<ModelType, true, false>;
    using Super::Super;
    using typename Super::RMSErrorType;
    using InputType  = typename ModelType::DataType;
    using OutputType = typename ModelType::ParameterType;
    using FlagType   = FlagType_; // Iterations

    virtual FitReturnType fit(const std::vector<QI_ARRAY(InputType)> &inputs,
                              typename ModelType::FixedArray const &  fixed,
                              typename ModelType::VaryingArray &      outputs,
                              typename ModelType::CovarArray *        cov,
                              RMSErrorType &                          rmse,
                              std::vector<QI_ARRAY(InputType)> &      point_residuals,
                              FlagType &                              flag,
                              const int                               block) const = 0;
};

template <typename ModelType, typename FlagType_ = int>
struct IndexedFitFunction : FitFunctionBase<ModelType, false, true> {
    using Super = FitFunctionBase<ModelType, false, true>;
    using Super::Super;
    using typename Super::RMSErrorType;
    using InputType    = typename ModelType::DataType;
    using OutputType   = typename ModelType::ParameterType;
    using FlagType     = FlagType_; // Iterations
    using SequenceType = typename ModelType::SequenceType;

    virtual FitReturnType fit(const std::vector<QI_ARRAY(InputType)> &inputs,
                              typename ModelType::FixedArray const &  fixed,
                              typename ModelType::VaryingArray &      outputs,
                              typename ModelType::CovarArray *        cov,
                              RMSErrorType &                          rmse,
                              std::vector<QI_ARRAY(InputType)> &      point_residuals,
                              FlagType &                              flag,
                              const itk::Index<3> &                   index) const = 0;
};

} // End namespace QI
