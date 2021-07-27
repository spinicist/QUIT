#pragma once

#include "FitFunction.h"

namespace QI {

template <typename ModelType_, typename FlagType_ = int> struct ScaledAutoDiffFit {
    using ModelType           = ModelType_;
    using InputType           = typename ModelType::DataType;
    using OutputType          = typename ModelType::ParameterType;
    using FlagType            = FlagType_; // Iterations
    using RMSErrorType        = double;
    static const bool Blocked = false;
    static const bool Indexed = false;

    ModelType model;
    long      input_size(long const &i) const { return model.input_size(i); }

    FitReturnType fit(std::vector<QI_ARRAY(InputType)> const &inputs,
                      typename ModelType::FixedArray const &  fixed,
                      typename ModelType::VaryingArray &      varying,
                      typename ModelType::DerivedArray &      derived,
                      typename ModelType::CovarArray *        cov,
                      RMSErrorType &                          rmse,
                      std::vector<QI_ARRAY(InputType)> &      residuals,
                      FlagType &                              iterations) const {
        const double &scale = inputs[0].maxCoeff();
        if (scale < std::numeric_limits<double>::epsilon()) {
            varying = ModelType::VaryingArray::Zero();
            rmse    = 0;
            return {false, "Maximum data value was not positive"};
        }
        const Eigen::ArrayXd data = inputs[0] / scale;
        ceres::Problem       problem;
        using Cost      = ModelCost<ModelType>;
        using AutoCost  = ceres::AutoDiffCostFunction<Cost, ceres::DYNAMIC, ModelType::NV>;
        auto *cost      = new Cost{this->model, fixed, data};
        auto *auto_cost = new AutoCost(cost, this->model.sequence.size());
        auto *loss      = new ceres::HuberLoss(1.0);
        problem.AddResidualBlock(auto_cost, loss, varying.data());
        for (int i = 0; i < ModelType::NV; i++) {
            problem.SetParameterLowerBound(varying.data(), i, this->model.bounds_lo[i]);
            problem.SetParameterUpperBound(varying.data(), i, this->model.bounds_hi[i]);
        }
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = 100;
        options.function_tolerance  = 1e-6;
        options.gradient_tolerance  = 1e-7;
        options.parameter_tolerance = 1e-5;
        options.logging_type        = ceres::SILENT;
        varying << this->model.start;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            return {false, summary.FullReport()};
        }
        iterations               = summary.iterations.size();
        Eigen::ArrayXd const rs  = (data - this->model.signal(varying, fixed));
        double const         var = rs.square().sum();
        rmse                     = sqrt(var / data.rows()) * scale;
        if (residuals.size() > 0) {
            residuals[0] = rs * scale;
        }
        if (cov) {
            QI::GetModelCovariance<ModelType>(
                problem, varying, var / (data.rows() - ModelType::NV), cov);
        }
        this->model.derived(varying, fixed, derived);
        varying[0] = varying[0] * scale;
        return {true, ""};
    }
};

} // namespace QI