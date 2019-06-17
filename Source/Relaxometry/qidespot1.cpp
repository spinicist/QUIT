/*
 *  qidespot1.cpp - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "ceres/ceres.h"
#include <Eigen/Core>
#include <array>

#include "Args.h"
#include "FitFunction.h"
#include "ImageIO.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "OnePoolSignals.h"
#include "SPGRSequence.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

using DESPOT1 = QI::Model<2, 1, QI::SPGRSequence>;
template <> std::array<const std::string, 2> DESPOT1::varying_names{{"PD"s, "T1"s}};
template <> std::array<const std::string, 1> DESPOT1::fixed_names{{"B1"s}};
template <> const QI_ARRAYN(double, 1) DESPOT1::fixed_defaults{1.0};
template <>
template <typename Derived>
auto DESPOT1::signal(const Eigen::ArrayBase<Derived> &v, const QI_ARRAYN(double, NF) & f) const
    -> QI_ARRAY(typename Derived::Scalar) {
    return QI::SPGRSignal(v[0], v[1], f[0], sequence);
}

using DESPOT1Fit = QI::FitFunction<DESPOT1>;

struct DESPOT1LLS : DESPOT1Fit {
    using DESPOT1Fit::DESPOT1Fit;
    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          const Eigen::ArrayXd &             fixed,
                          QI_ARRAYN(OutputType, DESPOT1::NV) & outputs,
                          RMSErrorType &               residual,
                          std::vector<Eigen::ArrayXd> &residuals,
                          FlagType &                   iterations) const override {
        const Eigen::ArrayXd &data = inputs[0];
        const double &        B1   = fixed[0];
        Eigen::ArrayXd        flip = model.sequence.FA * B1;
        Eigen::VectorXd       Y    = data / flip.sin();
        Eigen::MatrixXd       X(Y.rows(), 2);
        X.col(0) = data / flip.tan();
        X.col(1).setOnes();
        Eigen::Vector2d b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        outputs << QI::Clamp(b[1] / (1. - b[0]), model.bounds_lo[0], model.bounds_hi[0]),
            QI::Clamp(-model.sequence.TR / log(b[0]), model.bounds_lo[1], model.bounds_hi[1]);
        const Eigen::ArrayXd temp_residuals = data - model.signal(outputs, fixed);
        if (residuals.size() > 0) { // Residuals will only be allocated if the user asked for them
            residuals[0] = temp_residuals;
        }
        residual   = sqrt(temp_residuals.square().sum() / temp_residuals.rows());
        iterations = 1;
        return {true, ""};
    }
};

struct DESPOT1WLLS : DESPOT1Fit {
    using DESPOT1Fit::DESPOT1Fit;
    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          const Eigen::ArrayXd &             fixed,
                          QI_ARRAYN(OutputType, DESPOT1::NV) & outputs,
                          RMSErrorType &               residual,
                          std::vector<Eigen::ArrayXd> &residuals,
                          FlagType &                   iterations) const override {
        const Eigen::ArrayXd &data = inputs[0];
        const double &        B1   = fixed[0];
        Eigen::ArrayXd        flip = model.sequence.FA * B1;
        Eigen::VectorXd       Y    = data / flip.sin();
        Eigen::MatrixXd       X(Y.rows(), 2);
        X.col(0) = data / flip.tan();
        X.col(1).setOnes();
        Eigen::Vector2d b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        Eigen::Array2d  out{b[1] / (1. - b[0]), -model.sequence.TR / log(b[0])};
        for (iterations = 0; iterations < max_iterations; iterations++) {
            Eigen::VectorXd W =
                (flip.sin() / (1. - (exp(-model.sequence.TR / out[1]) * flip.cos()))).square();
            b = (X.transpose() * W.asDiagonal() * X)
                    .partialPivLu()
                    .solve(X.transpose() * W.asDiagonal() * Y);
            Eigen::Array2d newOut{b[1] / (1. - b[0]), -model.sequence.TR / log(b[0])};
            if (newOut.isApprox(out))
                break;
            else
                out = newOut;
        }
        // std::cout << "PD " << out[0] << " T1 " << out[1] << std::endl;
        outputs << QI::Clamp(out[0], model.bounds_lo[0], model.bounds_hi[0]),
            QI::Clamp(out[1], model.bounds_lo[1], model.bounds_hi[1]);
        const Eigen::ArrayXd temp_residuals = data - model.signal(outputs, fixed.cast<double>());
        if (residuals.size() > 0) { // Residuals will only be allocated if the user asked for them
            residuals[0] = temp_residuals;
        }
        residual = sqrt(temp_residuals.square().sum() / temp_residuals.rows());
        return {true, ""};
    }
};

struct DESPOT1NLLS : DESPOT1Fit {
    DESPOT1NLLS(DESPOT1 &m) : DESPOT1Fit(m) {
        // Do not let these go negative for non-linear fitting
        model.bounds_lo[0] = 1e-6;
        model.bounds_lo[1] = 1e-6;
    }

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          const Eigen::ArrayXd &             fixed,
                          QI_ARRAYN(OutputType, DESPOT1::NV) & p,
                          RMSErrorType &               residual,
                          std::vector<Eigen::ArrayXd> &residuals,
                          FlagType &                   iterations) const override {
        const double &scale = inputs[0].maxCoeff();
        if (scale < std::numeric_limits<double>::epsilon()) {
            p << 0.0, 0.0;
            residual = 0;
            return {false, "Maximum data value was not positive"};
        }
        const Eigen::ArrayXd data = inputs[0] / scale;
        p << 10., 1.;
        ceres::Problem problem;
        using Cost      = QI::ModelCost<DESPOT1>;
        using AutoCost  = ceres::AutoDiffCostFunction<Cost, ceres::DYNAMIC, DESPOT1::NV>;
        auto *cost      = new Cost(model, fixed, data);
        auto *auto_cost = new AutoCost(cost, model.sequence.size());
        problem.AddResidualBlock(auto_cost, NULL, p.data());
        problem.SetParameterLowerBound(p.data(), 0, model.bounds_lo[0] / scale);
        problem.SetParameterUpperBound(p.data(), 0, model.bounds_hi[0] / scale);
        problem.SetParameterLowerBound(p.data(), 1, model.bounds_lo[1]);
        problem.SetParameterUpperBound(p.data(), 1, model.bounds_hi[1]);
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = max_iterations;
        options.function_tolerance  = 1e-5;
        options.gradient_tolerance  = 1e-6;
        options.parameter_tolerance = 1e-4;
        options.logging_type        = ceres::SILENT;
        ceres::Solve(options, &problem, &summary);
        p[0] = p[0] * scale;
        if (!summary.IsSolutionUsable()) {
            return {false, summary.FullReport()};
        }
        iterations = summary.iterations.size();
        residual   = summary.final_cost * scale;
        if (residuals.size() > 0) {
            std::vector<double> r_temp(data.size());
            problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
            for (size_t i = 0; i < r_temp.size(); i++)
                residuals[0][i] = r_temp[i] * scale;
        }
        return {true, ""};
    }
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser(
        "Calculates T1 maps from SPGR data\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> spgr_path(parser, "SPGR FILE", "Path to SPGR data");
    QI_COMMON_ARGS;
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    args::ValueFlag<char> algorithm(parser, "ALGO", "Choose algorithm (l/w/n)", {'a', "algo"}, 'l');
    args::ValueFlag<int>  its(
        parser, "ITERS", "Max iterations for WLLS/NLLS (default 15)", {'i', "its"}, 15);
    args::ValueFlag<float> clampPD(parser,
                                   "CLAMP PD",
                                   "Clamp PD between 0 and value",
                                   {'p', "clampPD"},
                                   std::numeric_limits<float>::infinity());
    args::ValueFlag<float> clampT1(parser,
                                   "CLAMP T1",
                                   "Clamp T1 between 0 and value",
                                   {'t', "clampT1"},
                                   std::numeric_limits<float>::infinity());
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::CheckPos(spgr_path);
    QI::Log(verbose, "Reading sequence information");
    json input        = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    auto spgrSequence = input.at("SPGR").get<QI::SPGRSequence>();

    DESPOT1 model{spgrSequence};
    if (simulate) {
        QI::SimulateModel<DESPOT1, false>(
            input, model, {B1.Get()}, {spgr_path.Get()}, verbose, simulate.Get());
    } else {
        DESPOT1Fit *d1 = nullptr;
        switch (algorithm.Get()) {
        case 'l':
            d1 = new DESPOT1LLS(model);
            QI::Log(verbose, "LLS algorithm selected.");
            break;
        case 'w':
            d1 = new DESPOT1WLLS(model);
            QI::Log(verbose, "WLLS algorithm selected.");
            break;
        case 'n':
            d1 = new DESPOT1NLLS(model);
            QI::Log(verbose, "NLLS algorithm selected.");
            break;
        default:
            QI::Fail("Unknown algorithm type: {}", algorithm.Get());
        }
        if (clampPD) {
            d1->model.bounds_lo[0] = 1e-6;
            d1->model.bounds_hi[0] = clampPD.Get();
        }
        if (clampT1) {
            d1->model.bounds_lo[1] = 1e-6;
            d1->model.bounds_hi[1] = clampT1.Get();
        }
        if (its)
            d1->max_iterations = its.Get();
        auto fit =
            QI::ModelFitFilter<DESPOT1Fit>::New(d1, verbose.Get(), resids.Get(), subregion.Get());
        fit->ReadInputs({QI::CheckPos(spgr_path)}, {B1.Get()}, mask.Get());
        fit->Update();
        fit->WriteOutputs(prefix.Get() + "D1_");
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
