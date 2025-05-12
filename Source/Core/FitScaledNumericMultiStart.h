#pragma once
// #define QI_DEBUG_BUILD 1

#include "FitFunction.h"

namespace QI {

template <typename ModelType, int NScale = 1, int MultiStartParameter = -1>
struct ScaledNumericDiffMultiStartFit : FitFunction<ModelType, int> {
    using Super = FitFunction<ModelType, int>;
    using Super::Super;
    using typename Super::RMSErrorType;
    using InputType              = typename ModelType::DataType;
    using OutputType             = typename ModelType::ParameterType;
    static int constexpr MSIndex = (ModelType::NV + MultiStartParameter) % ModelType::NV;

    ScaledNumericDiffMultiStartFit(ModelType &m, Eigen::ArrayXd const &s) : Super{m}, starts{s} {}

    Eigen::ArrayXd starts;

    // This has to match the function signature that will be called in ModelFitFilter (which depends
    // on Blocked/Indexed. The return type is a simple struct indicating success, and on failure
    // also the reason for failure
    QI::FitReturnType
    fit(std::vector<QI_ARRAY(InputType)> const &inputs,  // Input: signal data
        typename ModelType::FixedArray const   &fixed,   // Input: Fixed parameters
        typename ModelType::VaryingArray       &varying, // Output: Varying parameters
        typename ModelType::CovarArray         *cov,
        RMSErrorType                           &rmse,      // Output: root-mean-square error
        std::vector<QI_ARRAY(InputType)>       &residuals, // Optional output: point residuals
        int &iterations /* Usually iterations */) const override {
        // First scale down the raw data so that PD will be roughly the same magnitude as other
        // parameters This is important for numerical stability in the optimiser
        InputType scale = inputs[0].maxCoeff();
        if (scale < std::numeric_limits<InputType>::epsilon()) {
            varying = ModelType::VaryingArray::Zero();
            rmse    = 0.0;
            return {false, "Maximum data value was not positive"};
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
                              this->model.input_size(0));
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
        options.max_num_iterations  = 30;
        options.function_tolerance  = 1e-6;
        options.gradient_tolerance  = 1e-7;
        options.parameter_tolerance = 1e-5;
        options.logging_type        = ceres::SILENT;

        double                  best_cost = std::numeric_limits<double>::infinity();
        typename ModelType::VaryingArray best_v;
        for (auto const s : this->starts) {
            varying          = this->model.start;
            varying(MSIndex) = s;

            ceres::Solve(options, &problem, &summary);
#ifdef QI_DEBUG_BUILD
            fmt::print(stderr, "Summary\n{}\n", summary.FullReport());
#endif
            if (!summary.IsSolutionUsable()) {
                return {false, summary.FullReport()};
            }

            double c;

            problem.Evaluate(ceres::Problem::EvaluateOptions(), &c, NULL, NULL, NULL);
            if (c < best_cost) {
                best_v     = varying;
                best_cost  = c;
                iterations = summary.iterations.size();
            }
        }
        varying = best_v;
        InputType              var;
        std::vector<InputType> rs(data.size());
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

} // namespace QI