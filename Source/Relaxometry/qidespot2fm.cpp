/*
 *  despot2fm.cpp
 *
 *  Created by Tobias Wood on 2015/06/03.
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
#include "SSFPSequence.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

using FMModel = QI::Model<3, 2, QI::SSFPSequence>;
template <> std::array<const std::string, 3> FMModel::varying_names{{"PD"s, "T2"s, "f0"s}};
template <> std::array<const std::string, 2> FMModel::fixed_names{{"T1"s, "B1"s}};
template <> const QI_ARRAYN(double, 2) FMModel::fixed_defaults{1.0, 1.0};
template <>
template <typename Derived>
auto FMModel::signal(const Eigen::ArrayBase<Derived> &v, const QI_ARRAYN(double, NF) & f) const
    -> QI_ARRAY(typename Derived::Scalar) {
    using T                      = typename Derived::Scalar;
    const T &     PD             = v[0];
    const T &     T2             = v[1];
    const T &     f0             = v[2];
    const double &T1             = f[0];
    const double &B1             = f[1];
    const double  E1             = exp(-sequence.TR / T1);
    const T       E2             = exp(-sequence.TR / T2);
    const T       psi            = 2. * M_PI * f0 * sequence.TR;
    const QI_ARRAY(double) alpha = sequence.FA * B1;
    const QI_ARRAY(T) d          = (1. - E1 * E2 * E2 - (E1 - E2 * E2) * cos(alpha));
    const QI_ARRAY(T) G          = -PD * (1. - E1) * sin(alpha) / d;
    const QI_ARRAY(T) b          = E2 * (1. - E1) * (1. + cos(alpha)) / d;

    const QI_ARRAY(T) theta  = sequence.PhaseInc + psi;
    const QI_ARRAY(T) cos_th = cos(theta);
    const QI_ARRAY(T) sin_th = sin(theta);
    const T cos_psi          = cos(psi);
    const T sin_psi          = sin(psi);
    const QI_ARRAY(T) re_m =
        (cos_psi - E2 * (cos_th * cos_psi - sin_th * sin_psi)) * G / (1.0 - b * cos_th);
    const QI_ARRAY(T) im_m =
        (sin_psi - E2 * (cos_th * sin_psi + sin_th * cos_psi)) * G / (1.0 - b * cos_th);
    return sqrt(re_m.square() + im_m.square());
}

using FMFit = QI::FitFunction<FMModel>;

struct FMNLLS : FMFit {
    using FMFit::FMFit;
    bool              asymmetric = false;
    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs, const Eigen::ArrayXd &fixed,
                          QI_ARRAYN(OutputType, FMModel::NV) & bestP, ResidualType &       residual,
                          std::vector<Eigen::ArrayXd> &residuals,
                          FlagType &                   iterations) const override {
        const double &T1 = fixed[0];
        if (std::isfinite(T1) && (T1 > model.sequence.TR)) {
            // Improve scaling by dividing the PD down to something sensible.
            // This gets scaled back up at the end.
            const double         scale = inputs[0].maxCoeff();
            const Eigen::ArrayXd data  = inputs[0] / scale;

            std::vector<double> f0_starts = {0, 0.4 / model.sequence.TR};
            if (this->asymmetric) {
                f0_starts.push_back(0.2 / model.sequence.TR);
                f0_starts.push_back(-0.2 / model.sequence.TR);
                f0_starts.push_back(-0.4 / model.sequence.TR);
            }

            double         best = std::numeric_limits<double>::infinity();
            Eigen::Array3d p;
            ceres::Problem problem;
            using Cost      = QI::ModelCost<FMModel>;
            using AutoCost  = ceres::AutoDiffCostFunction<Cost, ceres::DYNAMIC, FMModel::NV>;
            auto *cost      = new Cost(model, fixed, data);
            auto *auto_cost = new AutoCost(cost, model.sequence.size());
            problem.AddResidualBlock(auto_cost, NULL, p.data());
            problem.SetParameterLowerBound(p.data(), 0, 1.);
            problem.SetParameterLowerBound(p.data(), 1, model.sequence.TR);
            problem.SetParameterUpperBound(p.data(), 1, T1);
            if (this->asymmetric) {
                problem.SetParameterLowerBound(p.data(), 2, -0.5 / model.sequence.TR);
            } else {
                problem.SetParameterLowerBound(p.data(), 2, 0.0);
            }
            problem.SetParameterUpperBound(p.data(), 2, 0.5 / model.sequence.TR);
            ceres::Solver::Options options;
            ceres::Solver::Summary summary;
            options.max_num_iterations  = max_iterations;
            options.function_tolerance  = 1e-6;
            options.gradient_tolerance  = 1e-7;
            options.parameter_tolerance = 1e-5;
            options.logging_type        = ceres::SILENT;
            for (const double &f0 : f0_starts) {
                p = {
                    5., std::max(0.1 * T1, 1.5 * model.sequence.TR),
                    f0}; // Yarnykh gives T2 = 0.045 * T1 in brain, but best to overestimate for CSF
                ceres::Solve(options, &problem, &summary);
                if (!summary.IsSolutionUsable()) {
                    return std::make_tuple(false, summary.FullReport());
                }
                double r = summary.final_cost;
                if (r < best) {
                    best       = r;
                    bestP      = p;
                    residual   = summary.final_cost * scale;
                    iterations = summary.iterations.size();
                }
            }
            bestP[0] = bestP[0] * scale;
            if (!summary.IsSolutionUsable()) {
                return std::make_tuple(false, summary.FullReport());
            }
            if (residuals.size() > 0) {
                std::vector<double> r_temp(data.size());
                p = bestP; // Make sure the correct parameters are in the block
                problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
                for (size_t i = 0; i < r_temp.size(); i++)
                    residuals[i] = r_temp[i] * scale;
            }
        } else {
            bestP << 0.0, 0.0, 0.0;
            residual   = 0;
            iterations = 0;
            return std::make_tuple(false, "T1 was either infinite or shorter than TR");
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
        "Calculates a T2 map from SSFP data and a T1 map.\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> t1_path(parser, "T1_MAP", "Input T1 map");
    args::Positional<std::string> ssfp_path(parser, "SSFP_FILE", "Input SSFP file");

    args::HelpFlag       help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames",
                                        {'o', "out"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask",
                                      {'m', "mask"});
    args::Flag asym(parser, "ASYM", "Fit +/- off-resonance frequency", {'A', "asym"});
    args::ValueFlag<std::string> subregion(
        parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK",
        {'s', "subregion"});
    args::Flag resids(parser, "RESIDS", "Write out residuals for each data-point", {'r', "resids"});
    args::ValueFlag<int> its(parser, "ITERS", "Max iterations for NLLS (default 75)", {'i', "its"},
                             75);
    args::ValueFlag<std::string> seq_arg(parser, "FILE",
                                         "Read JSON input from file instead of stdin", {"file"});
    args::ValueFlag<float>       simulate(
        parser, "SIMULATE", "Simulate sequence instead of fitting model (argument is noise level)",
        {"simulate"}, 0.0);
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI::Log(verbose, "Reading sequence information");
    rapidjson::Document input = seq_arg ? QI::ReadJSON(seq_arg.Get()) : QI::ReadJSON(std::cin);
    QI::SSFPSequence    ssfp(QI::GetMember(input, "SSFP"));
    FMModel             model{ssfp};
    if (simulate) {
        QI::SimulateModel<FMModel, false>(input, model, {QI::CheckPos(t1_path), B1.Get()},
                                          {QI::CheckPos(ssfp_path)}, verbose, simulate.Get());
    } else {
        FMNLLS fm{model};
        fm.max_iterations = its.Get();
        fm.asymmetric     = asym.Get();
        auto fit_filter   = QI::ModelFitFilter<FMNLLS>::New(&fm, verbose, resids);
        fit_filter->SetInput(0, QI::ReadVectorImage(QI::CheckPos(ssfp_path), verbose));
        fit_filter->SetFixed(0, QI::ReadImage(QI::CheckPos(t1_path), verbose));
        if (B1)
            fit_filter->SetFixed(1, QI::ReadImage(B1.Get(), verbose));
        if (mask)
            fit_filter->SetMask(QI::ReadImage(mask.Get(), verbose));
        if (subregion)
            fit_filter->SetSubregion(QI::RegionArg(args::get(subregion)));
        fit_filter->Update();
        std::string outPrefix = outarg.Get() + "FM_";
        for (int i = 0; i < model.NV; i++) {
            QI::WriteImage(fit_filter->GetOutput(i),
                           outPrefix + fm.model.varying_names.at(i) + QI::OutExt(), verbose);
        }
        QI::WriteImage(fit_filter->GetResidualOutput(), outPrefix + "residual" + QI::OutExt(),
                       verbose);
        if (resids) {
            QI::WriteVectorImage(fit_filter->GetResidualsOutput(0),
                                 outPrefix + "all_residuals" + QI::OutExt(), verbose);
        }
        if (its) {
            QI::WriteImage(fit_filter->GetFlagOutput(), outPrefix + "iterations" + QI::OutExt(),
                           verbose);
        }
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
