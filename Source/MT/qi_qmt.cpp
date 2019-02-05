/*
 *  qi_qmt.cpp - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2018 Tobias Wood, Samuel Hurley, Erika Raven
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "ceres/ceres.h"
#include <Eigen/Core>

// #define QI_DEBUG_BUILD 1
#include "Args.h"
#include "FitFunction.h"
#include "ImageIO.h"
#include "Lineshape.h"
#include "MTSatSequence.h"
#include "Macro.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

struct RamaniReducedModel {
    using SequenceType  = QI::MTSatSequence;
    using DataType      = double;
    using ParameterType = double;

    static constexpr int NV = 5;
    static constexpr int ND = 0;
    static constexpr int NF = 2;

    using VaryingArray = QI_ARRAYN(ParameterType, NV);
    using FixedArray   = QI_ARRAYN(ParameterType, NF);
    static std::array<const std::string, NV> varying_names;
    static std::array<const std::string, NF> fixed_names;
    static const FixedArray                  fixed_defaults;

    SequenceType &                             sequence;
    const QI::Lineshapes                       lineshape;
    const std::shared_ptr<QI::InterpLineshape> interp = nullptr;

    VaryingArray bounds_lo;
    VaryingArray bounds_hi;
    VaryingArray start;

    RamaniReducedModel(SequenceType &s, const QI::Lineshapes ls,
                       const std::shared_ptr<QI::InterpLineshape> &i)
        : sequence{s}, lineshape{ls}, interp{i} {
        //{"R*M0a"s, "f/(R_a*(1-f))"s, "T2_b"s, "T1_a/T2_a"s, "gM0_a"s}
        start << 10.0, 0.1, 10.e-6, 25., 1.;
        // bounds_lo.setConstant(1.0e-12);
        // bounds_hi.setConstant(std::numeric_limits<double>::infinity());
        bounds_lo << 2, 1.e-12, 1.e-6, 1, 0.5;
        bounds_hi << 1e2, 1, 100.e-6, 50., 1.5;
        // bounds_hi << 100, 1, 50.e-6, 75, 1.5;
    }

    template <typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &v, const FixedArray &f) const
        -> QI_ARRAY(typename Derived::Scalar) {
        // Use Ramani's notation
        const auto &Rb      = 1.0; // Fix
        const auto &RM0a    = v[0];
        const auto &fterm   = v[1];
        const auto &T2b     = v[2];
        const auto &T1a_T2a = v[3];
        const auto &gM0a    = v[4];
        const auto &f0      = f[0];
        const auto &B1      = f[1];

        QI_ARRAY(typename Derived::Scalar) lsv;
        switch (lineshape) {
        case QI::Lineshapes::Gaussian:
            lsv = QI::Gaussian((sequence.sat_f0 + f0), T2b);
            break;
        case QI::Lineshapes::Lorentzian:
            lsv = QI::Lorentzian((sequence.sat_f0 + f0), T2b);
            break;
        case QI::Lineshapes::SuperLorentzian:
            lsv = QI::SuperLorentzian((sequence.sat_f0 + f0), T2b);
            break;
        case QI::Lineshapes::Interpolated:
            lsv = (*interp)((sequence.sat_f0 + f0), T2b);
            break;
        }

        const auto w_cwpe = (B1 * sequence.sat_angle / sequence.pulse.p1) *
                            sqrt(sequence.pulse.p2 / (sequence.pulse.Trf * sequence.TR));
        const auto R_rfb = M_PI * (w_cwpe * w_cwpe) * lsv;

        const auto S = gM0a * (Rb * RM0a * fterm + R_rfb + Rb + RM0a) /
                       ((RM0a * fterm) * (Rb + R_rfb) +
                        (1.0 + pow(w_cwpe / (2 * M_PI * sequence.sat_f0), 2.0) * T1a_T2a) *
                            (R_rfb + Rb + RM0a));
        return S;
    }
};
std::array<const std::string, 5> RamaniReducedModel::varying_names{
    {"R*M0a"s, "f/(R_a*(1-f))"s, "T2_b"s, "T1_a/T2_a"s, "gM0_a"s}};
std::array<const std::string, 2> RamaniReducedModel::fixed_names{{"f0"s, "B1"s}};
const QI_ARRAYN(double, 2) RamaniReducedModel::fixed_defaults{0.0, 1.0};

struct RamaniFullModel {
    using SequenceType  = QI::MTSatSequence;
    using DataType      = double;
    using ParameterType = double;

    static constexpr int NV = 6;
    static constexpr int ND = 0;
    static constexpr int NF = 3;

    using VaryingArray = QI_ARRAYN(ParameterType, NV);
    using FixedArray   = QI_ARRAYN(ParameterType, NF);
    static std::array<const std::string, NV> varying_names;
    static std::array<const std::string, NF> fixed_names;
    static const FixedArray                  fixed_defaults;

    VaryingArray bounds_lo = VaryingArray::Constant(1.0e-12);
    VaryingArray bounds_hi = VaryingArray::Constant(std::numeric_limits<ParameterType>::infinity());

    RamaniReducedModel  reduced;
    const SequenceType &sequence;

    RamaniFullModel(SequenceType &s, const QI::Lineshapes ls,
                    const std::shared_ptr<QI::InterpLineshape> &interp = nullptr)
        : reduced{s, ls, interp}, sequence{s} {}

    size_t num_outputs() const { return 2; }
    int    output_size(int i) {
        if (i == 0) {
            return sequence.size();
        } else {
            return 1;
        }
    }

    auto signals(const QI_ARRAY(double) & varying, const FixedArray &fixed) const
        -> std::vector<QI_ARRAY(double)> {
        RamaniReducedModel::VaryingArray reduced_v;
        const auto &                     PD   = varying[0];
        const auto &                     T1_f = varying[1];
        const auto &                     T2_f = varying[2];
        const auto &                     T1_b = 1.0; // Fixed
        const auto &                     T2_b = varying[3];
        const auto &                     k_bf = varying[4];
        const auto &                     f_b  = varying[5];

        const auto &R1a     = 1 / T1_f;
        const auto &R1b     = 1 / T1_b;
        const auto &RM0a    = k_bf * (1.0 - f_b) / f_b;
        const auto &RM0b    = k_bf;
        const auto &fterm   = f_b * T1_f / (1.0 - f_b);
        const auto &T1a_T2a = T1_f / T2_f;
        const auto &gM0a    = PD;

        reduced_v[0] = RM0a;
        reduced_v[1] = fterm;
        reduced_v[2] = T2_b;
        reduced_v[3] = T1a_T2a;
        reduced_v[4] = gM0a;

        Eigen::ArrayXd T1obs(1);
        T1obs << 2 / (RM0b + R1a + RM0a + R1b -
                      sqrt(pow(RM0b + R1a - RM0a - R1b, 2) + 4 * RM0a * RM0b));
        return {reduced.signal(reduced_v, fixed.head(2)), T1obs};
    }
};
std::array<const std::string, 6> RamaniFullModel::varying_names{
    {"PD"s, "T1_f"s, "T2_f"s, "T2_b"s, "k_bf"s, "f_b"s}};
std::array<const std::string, 3> RamaniFullModel::fixed_names{{"f0"s, "B1"s, "T1_app"}};
const QI_ARRAYN(double, 3) RamaniFullModel::fixed_defaults{0.0, 1.0, 1.0};

struct RamaniFitFunction : QI::FitFunction<RamaniFullModel> {
    using FitFunction::FitFunction;

    QI::FitReturnType fit(const std::vector<QI_ARRAY(InputType)> &inputs,
                          const Eigen::ArrayXd &                  fixed,
                          QI_ARRAYN(OutputType, RamaniFullModel::NV) & p, ResidualType &residual,
                          std::vector<QI_ARRAY(ResidualType)> &residuals,
                          FlagType &                           iterations) const override {
        const double &scale = inputs[0].maxCoeff();
        p                   = RamaniFullModel::VaryingArray::Zero();
        residual            = 0;
        if (scale < std::numeric_limits<double>::epsilon()) {
            return std::make_tuple(false, "Maximum data value was not positive");
        }
        const auto T1_obs = fixed[2];
        if (!std::isfinite(T1_obs) || (T1_obs < 1.e-12)) {
            return std::make_tuple(false, "T1 Observed was not finite and positive");
        }

        const Eigen::ArrayXd             data    = inputs[0] / scale;
        RamaniReducedModel::VaryingArray inner_p = model.reduced.start;
        QI_DBVEC(inner_p);
        ceres::Problem problem;
        using Cost      = QI::ModelCost<RamaniReducedModel>;
        using AutoCost  = ceres::AutoDiffCostFunction<Cost, ceres::DYNAMIC, RamaniReducedModel::NV>;
        auto *cost      = new Cost(model.reduced, fixed.head(2), data);
        auto *auto_cost = new AutoCost(cost, this->model.sequence.size());
        problem.AddResidualBlock(auto_cost, NULL, inner_p.data());
        for (int i = 0; i < 5; i++) {
            problem.SetParameterLowerBound(inner_p.data(), i, this->model.reduced.bounds_lo[i]);
            problem.SetParameterUpperBound(inner_p.data(), i, this->model.reduced.bounds_hi[i]);
        }
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = this->max_iterations;
        options.function_tolerance  = 1e-5;
        options.gradient_tolerance  = 1e-6;
        options.parameter_tolerance = 1e-4;
        options.logging_type        = ceres::SILENT;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            return std::make_tuple(false, summary.FullReport());
        }
        // Convert from the fitted parameters to useful ones
        const auto  R_obs   = 1 / T1_obs;
        const auto &Rb      = 1.0; // Fix
        const auto &RM0a    = inner_p[0];
        const auto &fterm   = inner_p[1];
        const auto &T2b     = inner_p[2];
        const auto &T1a_T2a = inner_p[3];
        const auto &gM0a    = inner_p[4];
        const auto  Ra      = ((Rb < R_obs) && (Rb - R_obs + RM0a) > 0)
                            ? R_obs
                            : R_obs / (1.0 + ((RM0a * fterm * (Rb - R_obs)) / (Rb - R_obs + RM0a)));
        const auto f    = fterm * Ra / (1.0 + fterm * Ra);
        const auto k_bf = RM0a * f / (1.0 - f);
        //{"R*M0a"s, "f/(R_a*(1-f))"s, "T2_b"s, "T1_a/T2_a"s, "gM0_a"s}
        //{"PD"s, "T1_f"s, "T2_f"s, "T2_b"s, "k_bf"s, "f_b"s}
        p[0] = gM0a * scale;
        p[1] = 1.0 / Ra;
        p[2] = p[1] / T1a_T2a;
        p[3] = T2b;
        p[4] = k_bf;
        p[5] = f;
        QI_DBVEC(inner_p);
        QI_DBVEC(p);
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
        "Calculates qMT maps from Gradient Echo Saturation data\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> mtsat_path(parser, "MTSAT FILE", "Path to MT-Sat data");
    args::Positional<std::string> T1(parser, "T1", "T1 map (seconds) file");
    args::HelpFlag                help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames",
                                        {'o', "out"});
    args::ValueFlag<std::string> f0(parser, "f0", "f0 map (Hz) file", {'f', "f0"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    args::ValueFlag<std::string> lineshape_arg(
        parser, "LINESHAPE",
        "Either Gaussian, Lorentzian, Superlorentzian, or a .json file generated by qi_lineshape",
        {'l', "lineshape"}, "Gaussian");
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask",
                                      {'m', "mask"});
    args::ValueFlag<std::string> subregion(
        parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK",
        {'s', "subregion"});
    args::Flag resids(parser, "RESIDS", "Write out residuals for each data-point", {'r', "resids"});
    args::ValueFlag<char> algorithm(parser, "ALGO", "Choose algorithm (l/w/n)", {'a', "algo"}, 'l');
    args::ValueFlag<std::string> seq_arg(parser, "FILE",
                                         "Read JSON input from file instead of stdin", {"file"});
    args::ValueFlag<float>       simulate(
        parser, "SIMULATE",
        "Simulate sequence instead of fit_filterting model (argument is noise level)", {"simulate"},
        0.0);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::CheckPos(mtsat_path);
    QI::Log(verbose, "Reading sequence information");
    rapidjson::Document input = seq_arg ? QI::ReadJSON(seq_arg.Get()) : QI::ReadJSON(std::cin);
    QI::MTSatSequence   mtsat_sequence(QI::GetMember(input, "MTSat"));
    QI::Lineshapes      lineshape;
    std::shared_ptr<QI::InterpLineshape> interp = nullptr;
    if (lineshape_arg.Get() == "Gaussian") {
        QI::Log(verbose, "Using a Gaussian lineshape");
        lineshape = QI::Lineshapes::Gaussian;
    } else if (lineshape_arg.Get() == "Lorentzian") {
        QI::Log(verbose, "Using a Lorentzian lineshape");
        lineshape = QI::Lineshapes::Lorentzian;
    } else if (lineshape_arg.Get() == "Superlorentzian") {
        QI::Log(verbose, "Using a Super-Lorentzian lineshape");
        lineshape = QI::Lineshapes::SuperLorentzian;
    } else {
        QI::Log(verbose, "Reading lineshape file: {}", lineshape_arg.Get());
        rapidjson::Document ls_file = QI::ReadJSON(lineshape_arg.Get());
        interp    = std::make_shared<QI::InterpLineshape>(QI::GetMember(ls_file, "lineshape"));
        lineshape = QI::Lineshapes::Interpolated;
    }

    RamaniFullModel model{mtsat_sequence, lineshape, interp};
    if (simulate) {
        QI::SimulateModel<RamaniFullModel, true>(input, model, {f0.Get(), B1.Get(), ""},
                                                 {mtsat_path.Get(), T1.Get()}, verbose,
                                                 simulate.Get());
    } else {
        RamaniFitFunction fit{model};

        auto fit_filter = QI::ModelFitFilter<RamaniFitFunction>::New(&fit, verbose, resids);
        fit_filter->SetInput(0, QI::ReadVectorImage(mtsat_path.Get(), verbose));
        if (f0)
            fit_filter->SetFixed(0, QI::ReadImage(f0.Get(), verbose));
        if (B1)
            fit_filter->SetFixed(1, QI::ReadImage(B1.Get(), verbose));
        if (T1)
            fit_filter->SetFixed(2, QI::ReadImage(T1.Get(), verbose));
        if (mask)
            fit_filter->SetMask(QI::ReadImage(mask.Get(), verbose));
        if (subregion)
            fit_filter->SetSubregion(QI::RegionArg(args::get(subregion)));
        fit_filter->Update();
        std::string outPrefix = outarg.Get() + "QMT_";
        for (int i = 0; i < RamaniFullModel::NV; i++) {
            QI::WriteImage(fit_filter->GetOutput(i),
                           outPrefix + RamaniFullModel::varying_names.at(i) + QI::OutExt());
        }
        QI::WriteImage(fit_filter->GetResidualOutput(), outPrefix + "residual" + QI::OutExt());
        if (resids) {
            QI::WriteVectorImage(fit_filter->GetResidualsOutput(0),
                                 outPrefix + "all_residuals" + QI::OutExt());
        }
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
