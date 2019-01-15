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
#include <iostream>

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

using MultiEcho = QI::Model<2, 0, QI::MultiEchoBase>;
template <>
template <typename Derived>
auto MultiEcho::signal(const Eigen::ArrayBase<Derived> &p, const QI_ARRAYN(double, 0) &
                       /*Unused*/) const -> QI_ARRAY(typename Derived::Scalar) {
    using T     = typename Derived::Scalar;
    const T &PD = p[0];
    const T &T2 = p[1];
    return PD * exp(-sequence.TE / T2);
}
template <> std::array<const std::string, 2> MultiEcho::varying_names{{"PD"s, "T2"s}};
template <> std::array<const std::string, 0> MultiEcho::fixed_names{{}};
template <> const QI_ARRAYN(double, 0) MultiEcho::fixed_defaults{};

using MultiEchoFit = QI::BlockFitFunction<MultiEcho>;

struct MultiEchoLogLin : MultiEchoFit {
    using MultiEchoFit::MultiEchoFit;
    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs, const Eigen::ArrayXd &fixed,
                          QI_ARRAYN(OutputType, MultiEcho::NV) & outputs, ResidualType &   residual,
                          std::vector<Eigen::ArrayXd> &residuals, FlagType &iterations,
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
        return std::make_tuple(true, "");
    }
};

struct MultiEchoARLO : MultiEchoFit {
    using MultiEchoFit::MultiEchoFit;
    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs, const Eigen::ArrayXd &fixed,
                          QI_ARRAYN(OutputType, MultiEcho::NV) & outputs, ResidualType &   residual,
                          std::vector<Eigen::ArrayXd> &residuals, FlagType &iterations,
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
        return std::make_tuple(true, "");
    }
};

struct MultiEchoNLLS : MultiEchoFit {
    using MultiEchoFit::MultiEchoFit;
    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs, const Eigen::ArrayXd &fixed,
                          QI_ARRAYN(OutputType, MultiEcho::NV) & p, ResidualType &         residual,
                          std::vector<Eigen::ArrayXd> &residuals, FlagType &iterations,
                          const int /*Unused*/) const override {
        const double &scale = inputs[0].maxCoeff();
        if (scale < std::numeric_limits<double>::epsilon()) {
            p << 0.0, 0.0;
            residual = 0;
            return std::make_tuple(false, "Maximum data value was not positive");
        }
        const Eigen::ArrayXd data = inputs[0] / scale;
        p << 10., 0.05;
        ceres::Problem problem;
        using Cost      = QI::ModelCost<MultiEcho>;
        using AutoCost  = ceres::AutoDiffCostFunction<Cost, ceres::DYNAMIC, MultiEcho::NV>;
        auto *cost      = new Cost(model, fixed, data);
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
        p[0] = p[0] * scale;
        if (!summary.IsSolutionUsable()) {
            return std::make_tuple(false, summary.FullReport());
        }
        iterations = summary.iterations.size();
        residual   = summary.final_cost * scale;
        if (residuals.size() > 0) {
            std::vector<double> r_temp(data.size());
            problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
            for (size_t i = 0; i < r_temp.size(); i++)
                residuals[0][i] = r_temp[i] * scale;
        }
        return std::make_tuple(true, "");
    }
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser(
        "Calculates T2/T2* maps from multi-echo data\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "INPUT FILE", "Input multi-echo data");
    args::HelpFlag                help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames",
                                        {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask",
                                      {'m', "mask"});
    args::ValueFlag<std::string> subregion(
        parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK",
        {'s', "subregion"});
    args::Flag resids(parser, "RESIDS", "Write out residuals for each data-point", {'r', "resids"});
    args::ValueFlag<char> algorithm(parser, "ALGO", "Choose algorithm (l/a/n)", {'a', "algo"}, 'l');
    args::ValueFlag<std::string> file(parser, "FILE", "Read JSON input from file instead of stdin",
                                      {"file"});
    args::ValueFlag<float>       simulate(
        parser, "SIMULATE", "Simulate sequence instead of fitting model (argument is noise level)",
        {"simulate"}, 0.0);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::CheckPos(input_path);
    QI_LOG(verbose, "Reading sequence parameters");
    rapidjson::Document input = file ? QI::ReadJSON(file.Get()) : QI::ReadJSON(std::cin);
    QI::MultiEchoBase * sequence;
    if (input.HasMember("MultiEcho")) {
        sequence = new QI::MultiEchoSequence(input["MultiEcho"]);
    } else if (input.HasMember("MultiEchoFlex")) {
        sequence = new QI::MultiEchoFlexSequence(input["MultiEchoFlex"]);
    } else {
        QI_FAIL("Could not find MultiEcho or MultiEchoFlex in JSON input");
    }
    MultiEcho model{*sequence};
    if (simulate) {
        QI::SimulateModel<MultiEcho, false>(input, model, {}, {QI::CheckPos(input_path)}, verbose,
                                            simulate.Get());
    } else {
        MultiEchoFit *me = nullptr;
        switch (algorithm.Get()) {
        case 'l':
            me = new MultiEchoLogLin(model);
            QI_LOG(verbose, "LogLin algorithm selected.");
            break;
        case 'a':
            me = new MultiEchoARLO(model);
            QI_LOG(verbose, "ARLO algorithm selected.");
            break;
        case 'n':
            me = new MultiEchoNLLS(model);
            QI_LOG(verbose, "Non-linear algorithm (Levenberg Marquardt) selected.");
            break;
        default:
            QI_FAIL("Unknown algorithm type " << algorithm.Get());
        }
        auto fit   = itk::ModelFitFilter<MultiEchoFit>::New(me, verbose, resids);
        auto input = QI::ReadVectorImage(input_path.Get(), verbose);
        fit->SetInput(0, input);
        const int nvols = input->GetNumberOfComponentsPerPixel();
        if (nvols % sequence->size() == 0) {
            const int nblocks = nvols / sequence->size();
            fit->SetBlocks(nblocks);
        } else {
            QI_FAIL("Input size is not a multiple of the sequence size");
        }
        if (mask)
            fit->SetMask(QI::ReadImage(mask.Get(), verbose));
        if (subregion)
            fit->SetSubregion(QI::RegionArg(subregion.Get()));
        fit->Update();
        std::string outPrefix = outarg.Get() + "ME_";
        for (int i = 0; i < model.NV; i++) {
            QI::WriteVectorImage(fit->GetOutput(i),
                                 outPrefix + me->model.varying_names.at(i) + QI::OutExt());
        }
        QI::WriteVectorImage(fit->GetResidualOutput(), outPrefix + "residual" + QI::OutExt());
        if (resids) {
            QI::WriteVectorImage(fit->GetResidualsOutput(0),
                                 outPrefix + "all_residuals" + QI::OutExt());
        }
        QI_LOG(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
