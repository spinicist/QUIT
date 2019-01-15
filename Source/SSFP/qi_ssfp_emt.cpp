/*
 *  qi_esmt.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "ceres/ceres.h"
#include <Eigen/Core>
#include <array>
#include <iostream>

#include "Args.h"
#include "ImageIO.h"
#include "ModelFitFilter.h"
#include "SSFPSequence.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

struct EMTModel {
    using DataType      = double;
    using ParameterType = double;
    using SequenceType  = QI::SSFPMTSequence;

    static constexpr int NV = 5;
    static constexpr int ND = 0;
    static constexpr int NF = 2;

    static std::array<const std::string, NV> varying_names;
    static std::array<const std::string, NF> fixed_names;
    static const QI_ARRAYN(double, NF) fixed_defaults;
    const SequenceType &sequence;
    double              T2_b = 12.e-6;

    size_t num_outputs() const { return 3; }
    int    output_size(int /* Unused */) { return sequence.size(); }

    template <typename Derived>
    auto signals(const Eigen::ArrayBase<Derived> &v, const QI_ARRAYN(double, NF) & f) const
        -> std::vector<QI_ARRAY(typename Derived::Scalar)> {
        using T             = typename Derived::Scalar;
        using ArrayXT       = Eigen::Array<T, Eigen::Dynamic, 1>;
        const T &     M0    = v[0];
        const T &     F     = v[1] / (1.0 - v[1]); // Convert from f_b to F
        const T &     k_bf  = v[2];
        const T &     T1_f  = v[3];
        const T &     T2_f  = v[4];
        const T &     T1_b  = T1_f;
        const double &f0_Hz = f[0];
        const double &B1    = f[1];

        const ArrayXT E1f   = (-sequence.TR / T1_f).exp();
        const ArrayXT E2_f  = (-sequence.TR / T2_f).exp();
        const ArrayXT E2_fe = (-sequence.TR / (2.0 * T2_f)).exp();
        const T       k_fb  = (F > 0.0) ? (k_bf / F) : T(0.0);
        const ArrayXT E1_b  = (-sequence.TR / T1_b).exp();
        const ArrayXT fk    = (-sequence.TR * (k_bf + k_fb)).exp();

        const double G_gauss =
            (T2_b / sqrt(2. * M_PI)) * exp(-pow(2. * M_PI * f0_Hz * T2_b, 2) / 2.0);
        const Eigen::ArrayXd WT =
            M_PI * (B1 * B1 * sequence.intB1) *
            G_gauss; // # Product of W and Trf to save a division and multiplication
        const Eigen::ArrayXd fw = (-WT).exp();
        const ArrayXT        A  = 1.0 + F - fw * E1_b * (F + fk);
        const ArrayXT        B  = 1.0 + fk * (F - fw * E1_b * (F + 1.0));
        const ArrayXT        C  = F * (1.0 - E1_b) * (1.0 - fk);

        const ArrayXT denom = (A - B * E1f * cos(B1 * sequence.FA) -
                               (E2_f * E2_f) * (B * E1f - A * cos(B1 * sequence.FA)));
        const ArrayXT Gp    = M0 * E2_fe * (sin(B1 * sequence.FA) * ((1.0 - E1f) * B + C)) / denom;
        const ArrayXT bp    = (E2_f * (A - B * E1f) * (1.0 + cos(B1 * sequence.FA))) / denom;
        const ArrayXT ap    = E2_f;

        return {Gp, ap, bp};
    }
};
std::array<const std::string, 5> EMTModel::varying_names{
    {"PD"s, "f_b"s, "k_bf"s, "T1_f"s, "T2_f"s}};
std::array<const std::string, 2> EMTModel::fixed_names{{"f0"s, "B1"s}};
const QI_ARRAYN(double, 2) EMTModel::fixed_defaults{0.0, 1.0};

struct EMTCost {
    const EMTModel &model;
    const QI_ARRAYN(double, EMTModel::NF) fixed;
    const QI_ARRAY(double) G, b;

    template <typename T> bool operator()(const T *const vin, T *rin) const {
        Eigen::Map<QI_ARRAY(T)>                            r(rin, G.rows() + b.rows());
        const Eigen::Map<const QI_ARRAYN(T, EMTModel::NV)> v(vin);
        const auto                                         signals = model.signals(v, fixed);
        r.head(G.rows())                                           = G - signals[0];
        r.tail(b.rows())                                           = b - signals[2];
        return true;
    }
};

struct EMTFit {
    static const bool Blocked = false;
    static const bool Indexed = false;
    using InputType           = double;
    using OutputType          = double;
    using ResidualType        = double;
    using FlagType            = int;
    using ModelType           = EMTModel;
    ModelType &model;

    int n_inputs() const { return 3; }
    int input_size(const int /* Unused */) const { return model.sequence.size(); }
    int n_fixed() const { return 2; }
    int n_outputs() const { return 5; }

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs, const Eigen::ArrayXd &fixed,
                          QI_ARRAYN(OutputType, EMTModel::NV) & p, ResidualType &          residual,
                          std::vector<Eigen::ArrayXd> &residuals, FlagType &iterations) const {
        const double          scale = inputs[0].mean();
        const Eigen::ArrayXd &G     = inputs[0] / scale;
        const Eigen::ArrayXd &a     = inputs[1];
        const Eigen::ArrayXd &b     = inputs[2];

        Eigen::ArrayXd T2_fs = (-model.sequence.TR / a.log());
        const double   T2_f  = T2_fs.mean(); // Different TRs so have to average afterwards

        auto *cost = new ceres::AutoDiffCostFunction<EMTCost, ceres::DYNAMIC, 5>(
            new EMTCost{model, fixed, G, b}, G.size() + b.size());
        ceres::LossFunction *loss = new ceres::HuberLoss(1.0);
        p << 15.0, 0.05, 5.0, 1.2, T2_f;
        ceres::Problem problem;
        problem.AddResidualBlock(cost, loss, p.data());
        problem.SetParameterLowerBound(p.data(), 0, 0.1);
        problem.SetParameterUpperBound(p.data(), 0, 20.0);
        problem.SetParameterLowerBound(p.data(), 1, 1e-6);
        problem.SetParameterUpperBound(p.data(), 1, 0.2 - 1e-6);
        problem.SetParameterLowerBound(p.data(), 2, 0.1);
        problem.SetParameterUpperBound(p.data(), 2, 20.0);
        problem.SetParameterLowerBound(p.data(), 3, 0.05);
        problem.SetParameterUpperBound(p.data(), 3, 5.0);
        problem.SetParameterLowerBound(p.data(), 4, T2_f * 0.999);
        problem.SetParameterUpperBound(p.data(), 4, T2_f * 1.001);
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = 100;
        options.function_tolerance  = 1e-7;
        options.gradient_tolerance  = 1e-8;
        options.parameter_tolerance = 1e-3;
        options.logging_type        = ceres::SILENT;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            return std::make_tuple(false, summary.FullReport());
        }
        p[0] *= scale;
        // f_b is converted to F internally so don't convert here
        residual = summary.final_cost;
        if (residuals.size() > 0) {
            std::vector<double> r_temp(G.size() + a.size() + b.size());
            problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
            for (int i = 0; i < G.size(); i++)
                residuals[i] = r_temp[i];
            Eigen::ArrayXd as = (-model.sequence.TR / T2_f).exp();
            for (int i = 0; i < a.size(); i++) {
                residuals[i + G.size()] = as[i] - a[i];
            }
            for (int i = 0; i < b.size(); i++) {
                residuals[i + G.size() + a.size()] = r_temp[i + G.size()];
            }
        }
        iterations = summary.iterations.size();
        return std::make_tuple(true, "");
    }
};

int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates qMT parameters from ellipse parameters.\nInputs are G, "
                                "a, b.\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> G_path(parser, "G_FILE", "Input G file");
    args::Positional<std::string> a_path(parser, "a_FILE", "Input a file");
    args::Positional<std::string> b_path(parser, "b_FILE", "Input b file");
    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag           debug(parser, "DEBUG", "Output debugging messages", {'d', "debug"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "PREFIX", "Add a prefix to output filenames",
                                        {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask",
                                      {'m', "mask"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio)", {'b', "B1"});
    args::ValueFlag<std::string> f0(parser, "f0", "f0 map (in Hertz)", {'f', "f0"});
    args::ValueFlag<double> T2_b_us(parser, "T2b", "T2 of bound pool (in microseconds, default 12)",
                                    {"T2b"}, 12);
    args::ValueFlag<std::string> subregion(
        parser, "REGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK",
        {'s', "subregion"});
    args::ValueFlag<std::string> seq_arg(parser, "FILE",
                                         "Read JSON input from file instead of stdin", {"file"});
    args::ValueFlag<float>       simulate(
        parser, "SIMULATE", "Simulate sequence instead of fitting model (argument is noise level)",
        {"simulate"}, 0.0);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::CheckPos(G_path);
    QI::CheckPos(a_path);
    QI::CheckPos(b_path);

    QI_LOG(verbose, "Reading sequence information");
    rapidjson::Document input = seq_arg ? QI::ReadJSON(seq_arg.Get()) : QI::ReadJSON(std::cin);
    QI::SSFPMTSequence  ssfp(QI::GetMember(input, "SSFPMT"));
    EMTModel            model{ssfp};
    if (simulate) {
        QI::SimulateModel<EMTModel, true>(input, model, {f0.Get(), B1.Get()},
                                          {G_path.Get(), a_path.Get(), b_path.Get()}, verbose,
                                          simulate.Get());
    } else {
        EMTFit fit{model};
        fit.model.T2_b  = T2_b_us.Get() * 1e-6;
        auto fit_filter = itk::ModelFitFilter<EMTFit>::New(&fit, verbose, false);
        fit_filter->SetInput(0, QI::ReadVectorImage(G_path.Get(), verbose));
        fit_filter->SetInput(1, QI::ReadVectorImage(a_path.Get(), verbose));
        fit_filter->SetInput(2, QI::ReadVectorImage(b_path.Get(), verbose));
        if (f0)
            fit_filter->SetFixed(0, QI::ReadImage(f0.Get(), verbose));
        if (B1)
            fit_filter->SetFixed(1, QI::ReadImage(B1.Get(), verbose));
        if (mask)
            fit_filter->SetMask(QI::ReadImage(mask.Get(), verbose));
        if (subregion)
            fit_filter->SetSubregion(QI::RegionArg(args::get(subregion)));
        fit_filter->Update();
        std::string outPrefix = outarg.Get() + "EMT_";
        for (int i = 0; i < EMTModel::NV; i++) {
            QI::WriteImage(fit_filter->GetOutput(i),
                           outPrefix + EMTModel::varying_names.at(i) + QI::OutExt());
        }
        QI_LOG(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
