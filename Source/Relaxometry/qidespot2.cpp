/*
 *  qidespot2.cpp
 *
 *  Created by Tobias Wood on 23/01/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
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
#include "SSFPSequence.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

struct DESPOT2 : QI::Model<2, 2, QI::SSFPSequence> {
    static std::array<const std::string, NV> varying_names;
    static std::array<const std::string, NF> fixed_names;
    static const QI_ARRAYN(double, NF) fixed_defaults;
    bool elliptical = false;

    template <typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &v, const QI_ARRAYN(double, NF) & f) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T          = typename Derived::Scalar;
        const T &     PD = v[0];
        const T &     T2 = v[1];
        const double &T1 = f[0];
        const double &B1 = f[1];
        const double  E1 = exp(-sequence.TR / T1);
        const T       E2 = exp(-sequence.TR / T2);

        const QI_ARRAY(double) alpha = sequence.FA * B1;
        const QI_ARRAY(T) denom = elliptical ? (1.0 - E1 * E2 * E2 - (E1 - E2 * E2) * cos(alpha))
                                             : (1.0 - E1 * E2 - (E1 - E2) * cos(alpha));
        const QI_ARRAY(T) numer = PD * sqrt(E2) * (1.0 - E1) * sin(alpha);
        return numer / denom;
    }
};
std::array<const std::string, 2> DESPOT2::varying_names{{"PD"s, "T2"s}};
std::array<const std::string, 2> DESPOT2::fixed_names{{"T1"s, "B1"s}};
const QI_ARRAYN(double, 2) DESPOT2::fixed_defaults{1.0, 1.0};

using DESPOT2Fit = QI::FitFunction<DESPOT2>;

struct DESPOT2LLS : DESPOT2Fit {
    using DESPOT2Fit::DESPOT2Fit;
    QI::FitReturnType fit(const std::vector<QI_ARRAY(InputType)> &inputs,
                          const Eigen::ArrayXd &                  fixed,
                          QI_ARRAYN(OutputType, DESPOT2::NV) & outputs,
                          ResidualType &                       residual,
                          std::vector<QI_ARRAY(ResidualType)> &residuals,
                          FlagType &                           iterations) const override {
        const Eigen::ArrayXd &data = inputs[0];
        const double &        T1   = fixed[0];
        const double &        B1   = fixed[1];
        const double &        TR   = model.sequence.TR;
        const double          E1   = exp(-TR / T1);
        double                PD, T2, E2;
        const Eigen::ArrayXd  angles = (model.sequence.FA * B1);

        Eigen::VectorXd Y = data / sin(angles);
        Eigen::MatrixXd X(Y.rows(), 2);
        X.col(0) = data / tan(angles);
        X.col(1).setOnes();
        Eigen::VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        if (model.elliptical) {
            T2 = 2. * TR / log((b[0] * E1 - 1.) / (b[0] - E1));
            E2 = exp(-TR / T2);
            PD = b[1] * (1. - E1 * E2 * E2) / (sqrt(E2) * (1. - E1));
        } else {
            T2 = TR / log((b[0] * E1 - 1.) / (b[0] - E1));
            E2 = exp(-TR / T2);
            PD = b[1] * (1. - E1 * E2) / (sqrt(E2) * (1. - E1));
        }
        outputs << QI::Clamp(PD, model.bounds_lo[0], model.bounds_hi[0]),
            QI::Clamp(T2, model.bounds_lo[1], model.bounds_hi[1]);
        Eigen::ArrayXd r = data - model.signal(outputs, fixed);
        if (residuals.size() > 0) { // Residuals will only be allocated if the user asked for them
            residuals[0] = r;
        }
        residual   = sqrt(r.square().sum() / r.rows());
        iterations = 1;
        return std::make_tuple(true, "");
    }
};

struct DESPOT2WLLS : DESPOT2Fit {
    using DESPOT2Fit::DESPOT2Fit;
    QI::FitReturnType fit(const std::vector<QI_ARRAY(InputType)> &inputs,
                          const Eigen::ArrayXd &                  fixed,
                          QI_ARRAYN(OutputType, DESPOT2::NV) & outputs,
                          ResidualType &                       residual,
                          std::vector<QI_ARRAY(ResidualType)> &residuals,
                          FlagType &                           iterations) const override {
        const Eigen::ArrayXd &data = inputs[0];
        const double &        T1   = fixed[0];
        const double &        B1   = fixed[1];
        const double          TR   = model.sequence.TR;
        const double          E1   = exp(-TR / T1);
        double                PD, T2, E2;
        const Eigen::ArrayXd  angles = (model.sequence.FA * B1);

        Eigen::VectorXd Y = data / angles.sin();
        Eigen::MatrixXd X(Y.rows(), 2);
        X.col(0) = data / angles.tan();
        X.col(1).setOnes();
        Eigen::VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        if (model.elliptical) {
            T2 = 2. * TR / log((b[0] * E1 - 1.) / (b[0] - E1));
            E2 = exp(-TR / T2);
            PD = b[1] * (1. - E1 * E2 * E2) / (sqrt(E2) * (1. - E1));
        } else {
            T2 = TR / log((b[0] * E1 - 1.) / (b[0] - E1));
            E2 = exp(-TR / T2);
            PD = b[1] * (1. - E1 * E2) / (1. - E1);
        }
        Eigen::VectorXd W(model.sequence.size());
        for (iterations = 0; iterations < max_iterations; iterations++) {
            if (model.elliptical) {
                W = ((1. - E1 * E2) * angles.sin() /
                     (1. - E1 * E2 * E2 - (E1 - E2 * E2) * angles.cos()))
                        .square();
            } else {
                W = ((1. - E1 * E2) * angles.sin() / (1. - E1 * E2 - (E1 - E2) * angles.cos()))
                        .square();
            }
            b = (X.transpose() * W.asDiagonal() * X)
                    .partialPivLu()
                    .solve(X.transpose() * W.asDiagonal() * Y);
            if (model.elliptical) {
                T2 = 2. * TR / log((b[0] * E1 - 1.) / (b[0] - E1));
                E2 = exp(-TR / T2);
                PD = b[1] * (1. - E1 * E2 * E2) / (sqrt(E2) * (1. - E1));
            } else {
                T2 = TR / log((b[0] * E1 - 1.) / (b[0] - E1));
                E2 = exp(-TR / T2);
                PD = b[1] * (1. - E1 * E2) / (1. - E1);
            }
        }
        outputs[0] = QI::Clamp(PD, model.bounds_lo[0], model.bounds_hi[0]);
        outputs[1] = QI::Clamp(T2, model.bounds_lo[1], model.bounds_hi[1]);
        Eigen::Array2d v{PD, T2};
        Eigen::ArrayXd r = data - model.signal(v, fixed);
        if (residuals.size() > 0) { // Residuals will only be allocated if the user asked for them
            residuals[0] = r;
        }
        residual = sqrt(r.square().sum() / r.rows());
        return std::make_tuple(true, "");
    }
};

struct DESPOT2NLLS : DESPOT2Fit {
    DESPOT2NLLS(DESPOT2 &m) : DESPOT2Fit{m} {
        model.bounds_lo[0] = 1e-6; // Don't go negative PD
        model.bounds_lo[1] = 1e-3; // The sqrt(E2) term goes crazy for T2 below this
    }

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          const Eigen::ArrayXd &             fixed,
                          QI_ARRAYN(OutputType, DESPOT2::NV) & p,
                          ResidualType &               residual,
                          std::vector<Eigen::ArrayXd> &residuals,
                          FlagType &                   iterations) const override {
        const double &scale = inputs[0].maxCoeff();
        if (scale < std::numeric_limits<double>::epsilon()) {
            p << 0.0, 0.0;
            residual = 0;
            return std::make_tuple(false, "Maximum data value was not positive");
        }
        const Eigen::ArrayXd data = inputs[0] / scale;
        p << 10., 0.1;
        ceres::Problem problem;
        using Cost      = QI::ModelCost<DESPOT2>;
        using AutoCost  = ceres::AutoDiffCostFunction<Cost, ceres::DYNAMIC, DESPOT2::NV>;
        auto *cost      = new Cost(model, fixed, data);
        auto *auto_cost = new AutoCost(cost, model.sequence.size());
        problem.AddResidualBlock(auto_cost, NULL, p.data());
        problem.SetParameterLowerBound(p.data(), 0, model.bounds_lo[0] / scale);
        problem.SetParameterUpperBound(p.data(), 0, model.bounds_hi[0] / scale);
        problem.SetParameterLowerBound(p.data(), 1, model.bounds_lo[1]);
        problem.SetParameterUpperBound(p.data(), 1,
                                       std::min(model.bounds_hi[1], fixed[0])); // T2 cannot be > T1
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
        "Calculates T2 maps from SSFP data\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> t1_path(parser, "T1 MAP", "Path to T1 map");
    args::Positional<std::string> ssfp_path(parser, "SSFP FILE", "Path to SSFP data");
    args::HelpFlag                help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames",
                                        {'o', "out"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask",
                                      {'m', "mask"});
    args::ValueFlag<std::string> subregion(
        parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK",
        {'s', "subregion"});
    args::Flag resids(parser, "RESIDS", "Write out residuals for each data-point", {'r', "resids"});
    args::ValueFlag<char> algorithm(parser, "ALGO", "Choose algorithm (l/w/n)", {'a', "algo"}, 'l');
    args::Flag           gs_arg(parser, "GS", "Data is band-free geometric solution / ellipse data",
                      {'g', "gs"});
    args::ValueFlag<int> its(parser, "ITERS", "Max iterations for WLLS/NLLS (default 15)",
                             {'i', "its"}, 15);
    args::ValueFlag<float>       clampPD(parser, "CLAMP PD", "Clamp PD between 0 and value",
                                   {'p', "clampPD"}, std::numeric_limits<float>::infinity());
    args::ValueFlag<float>       clampT2(parser, "CLAMP T2", "Clamp T2 between 0 and value",
                                   {'t', "clampT2"}, std::numeric_limits<float>::infinity());
    args::ValueFlag<std::string> seq_arg(parser, "FILE",
                                         "Read JSON input from file instead of stdin", {"file"});
    args::ValueFlag<float>       simulate(
        parser, "SIMULATE", "Simulate sequence instead of fitting model (argument is noise level)",
        {"simulate"}, 0.0);
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI::Log(verbose, "Reading sequence information");
    rapidjson::Document input = seq_arg ? QI::ReadJSON(seq_arg.Get()) : QI::ReadJSON(std::cin);
    QI::SSFPSequence    ssfp(QI::GetMember(input, "SSFP"));
    DESPOT2             model{ssfp};
    if (simulate) {
        if (gs_arg)
            model.elliptical = true;
        QI::SimulateModel<DESPOT2, false>(input, model, {QI::CheckPos(t1_path), B1.Get()},
                                          {QI::CheckPos(ssfp_path)}, verbose, simulate.Get());
    } else {
        DESPOT2Fit *d2 = nullptr;
        switch (algorithm.Get()) {
        case 'l':
            d2 = new DESPOT2LLS(model);
            QI::Log(verbose, "LLS algorithm selected.");
            break;
        case 'w':
            d2 = new DESPOT2WLLS(model);
            QI::Log(verbose, "WLLS algorithm selected.");
            break;
        case 'n':
            d2 = new DESPOT2NLLS(model);
            QI::Log(verbose, "NLLS algorithm selected.");
            break;
        }
        if (clampPD) {
            d2->model.bounds_lo[0] = 1e-6;
            d2->model.bounds_hi[0] = clampPD.Get();
        }
        if (clampT2) {
            d2->model.bounds_lo[1] = 1e-6;
            d2->model.bounds_hi[1] = clampT2.Get();
        }
        if (its)
            d2->max_iterations = its.Get();
        if (gs_arg) {
            QI::Log(verbose, "GS Mode selected");
            d2->model.elliptical = true;
        }
        auto fit = QI::ModelFitFilter<DESPOT2Fit>::New(d2, verbose, resids, subregion.Get());
        fit->ReadInputs({QI::CheckPos(ssfp_path)}, {QI::CheckPos(t1_path), B1.Get()}, mask.Get());
        fit->Update();
        fit->WriteOutputs(outarg.Get() + "D2_");
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
