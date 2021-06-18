/*
 *  multiecho.cpp
 *
 *  Created by Tobias Wood on 27/01/2015.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <array>

#include "ceres/ceres.h"
#include <Eigen/Core>

#include "Args.h"
#include "FitFunction.h"
#include "ImageIO.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "MultiEchoSequence.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

struct MultiEcho : QI::Model<double, double, 2, 0> {
    QI::MultiEchoSequence const &sequence;

    std::array<const std::string, 2> const varying_names{{"PD"s, "T2"s}};
    VaryingArray const                     start{10., 0.05};
    VaryingArray const                     bounds_lo{0.1, 0.001};
    VaryingArray const                     bounds_hi{100., 5.};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    template <typename Derived>
    auto signal(Eigen::ArrayBase<Derived> const &p, FixedArray const &
                /*Unused*/) const -> QI_ARRAY(typename Derived::Scalar) {
        using T     = typename Derived::Scalar;
        const T &PD = p[0];
        const T &T2 = p[1];
        return PD * exp(-sequence.TE / T2);
    }
};

using MultiEchoFit = QI::BlockFitFunction<MultiEcho>;

struct MultiEchoLogLin : MultiEchoFit {
    using MultiEchoFit::MultiEchoFit;
    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          MultiEcho::FixedArray const &      fixed,
                          MultiEcho::VaryingArray &          outputs,
                          MultiEcho::CovarArray * /* Unused */,
                          RMSErrorType &               residual,
                          std::vector<Eigen::ArrayXd> &residuals,
                          FlagType &                   iterations,
                          const int /*Unused*/) const override {
        const Eigen::ArrayXd &data = inputs[0];
        Eigen::MatrixXd       X(model.sequence.size(), 2);
        X.col(0) = model.sequence.TE;
        X.col(1).setOnes();
        Eigen::VectorXd Y = data.array().log();
        Eigen::VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        outputs << exp(b[1]), -1 / b[0];
        const Eigen::ArrayXd temp_residuals = data - model.signal(outputs, fixed);
        if (residuals.size() > 0) { // Residuals will only be allocated if the user asked for them
            residuals[0] = temp_residuals;
        }
        residual   = sqrt(temp_residuals.square().sum() / temp_residuals.rows());
        iterations = 1;
        return {true, ""};
    }
};

struct MultiEchoARLO : MultiEchoFit {
    using MultiEchoFit::MultiEchoFit;
    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          MultiEcho::FixedArray const &      fixed,
                          MultiEcho::VaryingArray &          outputs,
                          MultiEcho::CovarArray * /*Unused*/,
                          RMSErrorType &               residual,
                          std::vector<Eigen::ArrayXd> &residuals,
                          FlagType &                   iterations,
                          const int /*Unused*/) const override {
        const Eigen::ArrayXd &data   = inputs[0];
        const double          ESP    = model.sequence.TE[1] - model.sequence.TE[0];
        const double          dTE_3  = (ESP / 3);
        double                si2sum = 0, di2sum = 0, sidisum = 0;
        for (Eigen::Index i = 0; i < model.sequence.size() - 2; i++) {
            const double si = dTE_3 * (data(i) + 4 * data(i + 1) + data(i + 2));
            const double di = data(i) - data(i + 2);
            si2sum += si * si;
            di2sum += di * di;
            sidisum += si * di;
        }
        double T2 = (si2sum + dTE_3 * sidisum) / (dTE_3 * di2sum + sidisum);
        double PD = (data.array() / exp(-model.sequence.TE / T2)).mean();
        outputs << PD, T2;
        const Eigen::ArrayXd temp_residuals = data - model.signal(outputs, fixed);
        if (residuals.size() > 0) { // Residuals will only be allocated if the user asked for them
            residuals[0] = temp_residuals;
        }
        residual   = sqrt(temp_residuals.square().sum() / temp_residuals.rows());
        iterations = 1;
        return {true, ""};
    }
};

struct MultiEchoNLLS : MultiEchoFit {
    using MultiEchoFit::MultiEchoFit;
    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          MultiEcho::FixedArray const &      fixed,
                          MultiEcho::VaryingArray &          p,
                          MultiEcho::CovarArray *            cov,
                          RMSErrorType &                     rmse,
                          std::vector<Eigen::ArrayXd> &      residuals,
                          FlagType &                         iterations,
                          const int /*Unused*/) const override {
        const double &scale = inputs[0].maxCoeff();
        if (scale < std::numeric_limits<double>::epsilon()) {
            p << 0.0, 0.0;
            rmse = 0;
            return {false, "Maximum data value was not positive"};
        }
        const Eigen::ArrayXd data = inputs[0] / scale;
        p                         = model.start;
        ceres::Problem problem;
        using Cost      = QI::ModelCost<MultiEcho>;
        using AutoCost  = ceres::AutoDiffCostFunction<Cost, ceres::DYNAMIC, MultiEcho::NV>;
        auto *cost      = new Cost{model, fixed, data};
        auto *auto_cost = new AutoCost(cost, model.sequence.size());
        problem.AddResidualBlock(auto_cost, NULL, p.data());
        problem.SetParameterLowerBound(p.data(), 0, 1.0e-6);
        problem.SetParameterUpperBound(p.data(), 0, model.bounds_hi[0] / scale);
        problem.SetParameterLowerBound(p.data(), 1, 1.0e-3);
        problem.SetParameterUpperBound(p.data(), 1, model.bounds_hi[1]);
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = 50;
        options.function_tolerance  = 1e-5;
        options.gradient_tolerance  = 1e-6;
        options.parameter_tolerance = 1e-4;
        options.logging_type        = ceres::SILENT;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            return {false, summary.FullReport()};
        }
        iterations = summary.iterations.size();

        Eigen::ArrayXd const rs  = (data - model.signal(p, fixed));
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

//******************************************************************************
// Main
//******************************************************************************
int multiecho_main(args::Subparser &parser) {
    args::Positional<std::string> input_path(parser, "INPUT FILE", "Input multi-echo data");
    QI_COMMON_ARGS;
    args::ValueFlag<char> algorithm(parser, "ALGO", "Choose algorithm (l/a/n)", {'a', "algo"}, 'l');
    parser.Parse();
    QI::CheckPos(input_path);
    QI::Log(verbose, "Reading sequence parameters");
    json input    = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    auto sequence = input.at("MultiEcho").get<QI::MultiEchoSequence>();

    MultiEcho model{{}, sequence};
    if (simulate) {
        QI::SimulateModel<MultiEcho, false>(input,
                                            model,
                                            {},
                                            {QI::CheckPos(input_path)},
                                            mask.Get(),
                                            verbose,
                                            simulate.Get(),
                                            subregion.Get());
    } else {
        MultiEchoFit *me = nullptr;
        switch (algorithm.Get()) {
        case 'l':
            me = new MultiEchoLogLin(model);
            QI::Log(verbose, "LogLin algorithm selected.");
            break;
        case 'a':
            me = new MultiEchoARLO(model);
            QI::Log(verbose, "ARLO algorithm selected.");
            break;
        case 'n':
            me = new MultiEchoNLLS(model);
            QI::Log(verbose, "Non-linear algorithm (Levenberg Marquardt) selected.");
            break;
        default:
            QI::Fail("Unknown algorithm type {}", algorithm.Get());
        }
        auto fit =
            QI::ModelFitFilter<MultiEchoFit>::New(
                me, verbose, covar, resids, threads.Get(), subregion.Get());
        fit->ReadInputs({QI::CheckPos(input_path)}, {}, mask.Get());
        const int nvols = fit->GetInput(0)->GetNumberOfComponentsPerPixel();
        if (nvols % sequence.size() == 0) {
            const int nblocks = nvols / sequence.size();
            fit->SetBlocks(nblocks);
        } else {
            QI::Fail("Input size is not a multiple of the sequence size");
        }
        fit->Update();
        fit->WriteOutputs(prefix.Get() + "ME_");
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
