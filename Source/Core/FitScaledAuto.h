#pragma once

#include "FitFunction.h"

namespace QI {

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
        options.max_num_iterations  = 30;
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

} // namespace QI