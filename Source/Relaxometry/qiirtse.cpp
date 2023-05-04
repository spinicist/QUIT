/*
 *  qiirtse.cpp
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
#include "IRTSESequence.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

struct IRTSE : QI::Model<double, double, 2, 0> {
    QI::IRTSESequence const &sequence;

    std::array<const std::string, 2> const varying_names{{"PD"s, "T1"s}};
    VaryingArray const                     start{1., 0.1};
    VaryingArray const                     bounds_lo{0.01, 0.001};
    VaryingArray const                     bounds_hi{10000., 10.};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    template <typename Derived>
    auto signal(Eigen::ArrayBase<Derived> const &p, FixedArray const &
                /*Unused*/) const -> QI_ARRAY(typename Derived::Scalar) {
        using T     = typename Derived::Scalar;
        const T &PD = p[0];
        const T &T1 = p[1];
        auto M = PD * (1. - exp(-sequence.TI/T1)) + sequence.Q * PD*exp(-sequence.TI/T1)*((1. - exp(-sequence.TD1/T1))*exp(-sequence.TD2/T1)*cos(sequence.theta) + (1. - exp(-sequence.TD2/T1)));
        
        return M;
    }
};

using IRTSEFit = QI::BlockFitFunction<IRTSE>;

struct IRTSENLLS : IRTSEFit {
    using IRTSEFit::IRTSEFit;
    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          IRTSE::FixedArray const &      fixed,
                          IRTSE::VaryingArray &          p,
                          IRTSE::CovarArray *            cov,
                          RMSErrorType &                     rmse,
                          std::vector<Eigen::ArrayXd> &      residuals,
                          FlagType &                         iterations,
                          const int /*Unused*/) const override {
        
        const double scale = std::abs(inputs[0].maxCoeff());
        const Eigen::ArrayXd data = inputs[0] / scale;
        
        p                         = model.start;
        ceres::Problem problem;
        using Cost      = QI::ModelCost<IRTSE>;
        using AutoCost  = ceres::AutoDiffCostFunction<Cost, ceres::DYNAMIC, IRTSE::NV>;
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
int irtse_main(args::Subparser &parser) {
    args::Positional<std::string> input_path(parser, "INPUT FILE", "Input multi-TI data");
    QI_COMMON_ARGS;
    parser.Parse();
    QI::CheckPos(input_path);
    QI::Log(verbose, "Reading sequence parameters");
    json input    = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    auto sequence = input.at("IRTSE").get<QI::IRTSESequence>();
    IRTSE model{{}, sequence};
    if (simulate) {
        QI::SimulateModel<IRTSE, false>(input,
                                            model,
                                            {},
                                            {QI::CheckPos(input_path)},
                                            mask.Get(),
                                            verbose,
                                            simulate.Get(),
                                            threads.Get(),
                                            subregion.Get());
    } else {
        IRTSEFit *me = nullptr;

        me = new IRTSENLLS(model);
        QI::Log(verbose, "Non-linear algorithm (Levenberg Marquardt) selected.");
        
        auto fit =
            QI::ModelFitFilter<IRTSEFit>::New(
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
        fit->WriteOutputs(prefix.Get() + "IRTSE_");
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
