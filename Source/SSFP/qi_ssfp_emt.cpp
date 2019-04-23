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

// #define QI_DEBUG_BUILD 1

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

    static constexpr int NV = 4;
    static constexpr int ND = 0;
    static constexpr int NF = 3;

    SequenceType const &  sequence;
    Eigen::ArrayXd const &Wtot;
    double const &        T2_b = 12.e-6;

    std::array<std::string, NV> const varying_names{"PD"s, "f_b"s, "k_bf"s, "T1_f"s};
    std::array<std::string, NF> const fixed_names{"f0"s, "B1"s, "T2_f"s};
    QI_ARRAYN(double, NF) const fixed_defaults{0.0, 1.0, 0.1};

    size_t num_outputs() const { return 3; }
    int    output_size(int /* Unused */) { return sequence.size(); }

    template <typename Derived>
    auto signals(const Eigen::ArrayBase<Derived> &v, const QI_ARRAYN(double, NF) & f) const
        -> std::vector<QI_ARRAY(typename Derived::Scalar)> {
        using T             = typename Derived::Scalar;
        using ArrayXT       = Eigen::Array<T, Eigen::Dynamic, 1>;
        const T &     M0    = v[0];
        const T &     f_b   = v[1];
        const T &     f_f   = 1.0 - f_b;
        const T &     k_bf  = v[2];
        const T &     T1_f  = v[3];
        const T &     T1_b  = T1_f;
        double const &f0_Hz = f[0];
        double const &B1    = f[1];
        double const &T2_f  = f[2];

        const ArrayXT        E1f   = (-sequence.TR / T1_f).exp();
        const Eigen::ArrayXd E2_f  = (-sequence.TR / T2_f).exp();
        const Eigen::ArrayXd E2_fe = (-sequence.TR / (2.0 * T2_f)).exp();
        const T              k_fb  = (f_b > 0.0) ? (k_bf * f_f / f_b) : T(0.0);
        const ArrayXT        E1_b  = (-sequence.TR / T1_b).exp();
        const ArrayXT        Ek    = (-sequence.TR * (k_bf + k_fb)).exp();

        const double G_gauss =
            sqrt(1.0 / (2.0 * M_PI)) * T2_b * exp(-pow(2.0 * M_PI * f0_Hz * T2_b, 2.0) / 2.0);
        const Eigen::ArrayXd Ew = (-Wtot * B1 * B1 * M_PI * M_PI * G_gauss).exp();
        const ArrayXT        A  = 1.0 - Ew * E1_b * (f_b + f_f * Ek);
        const ArrayXT        B  = f_f - Ek * (Ew * E1_b - f_b);
        const ArrayXT        C  = f_b * (1.0 - E1_b) * (1.0 - Ek);

        const ArrayXT denom = A - B * E1f * cos(B1 * sequence.FA) -
                              (E2_f * E2_f) * (B * E1f - A * cos(B1 * sequence.FA));
        const ArrayXT G = M0 * E2_fe * (sin(B1 * sequence.FA) * (B * (1.0 - E1f) + C)) / denom;
        const ArrayXT b = (E2_f * (A - B * E1f) * (1.0 + cos(B1 * sequence.FA))) / denom;

        // Annoying hack for simulating data
        ArrayXT a(sequence.size());
        for (Eigen::Index i = 0; i < sequence.size(); i++) {
            a[i] = T(exp(-sequence.TR[i] / T2_f));
        }
        return {G, a, b};
    }
};

struct EMTCost {
    const EMTModel &model;
    const QI_ARRAYN(double, EMTModel::NF) fixed;
    const QI_ARRAY(double) G, b;

    template <typename T> bool operator()(const T *const vin, T *rin) const {
        Eigen::Map<QI_ARRAY(T)>                            r(rin, G.rows() + b.rows());
        const Eigen::Map<const QI_ARRAYN(T, EMTModel::NV)> v(vin);

        const auto signals = model.signals(v, fixed);
        r.head(G.rows())   = G - signals[0];
        r.tail(b.rows())   = b - signals[2];

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

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          const Eigen::ArrayXd &             fixed,
                          QI_ARRAYN(OutputType, EMTModel::NV) & p,
                          ResidualType &               residual,
                          std::vector<Eigen::ArrayXd> &residuals,
                          FlagType &                   iterations) const {
        const double          scale = inputs[0].mean();
        const Eigen::ArrayXd &G     = inputs[0] / scale;
        const Eigen::ArrayXd &a     = inputs[1];
        const Eigen::ArrayXd &b     = inputs[2];

        auto *cost = new ceres::AutoDiffCostFunction<EMTCost, ceres::DYNAMIC, 5>(
            new EMTCost{model, fixed, G, b}, G.size() + b.size());
        ceres::LossFunction *loss = new ceres::HuberLoss(1.0);
        p << 10.0, 0.05, 2.0, 1.0;
        ceres::Problem problem;
        problem.AddResidualBlock(cost, loss, p.data());
        problem.SetParameterLowerBound(p.data(), 0, 0.1);
        problem.SetParameterUpperBound(p.data(), 0, 20.0);
        problem.SetParameterLowerBound(p.data(), 1, 1e-6);
        problem.SetParameterUpperBound(p.data(), 1, 0.2 - 1e-6);
        problem.SetParameterLowerBound(p.data(), 2, 0.1);
        problem.SetParameterUpperBound(p.data(), 2, 10.0);
        problem.SetParameterLowerBound(p.data(), 3, 0.05);
        problem.SetParameterUpperBound(p.data(), 3, 5.0);
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = 100;
        options.function_tolerance  = 1e-7;
        options.gradient_tolerance  = 1e-8;
        options.parameter_tolerance = 1e-3;
        options.logging_type        = ceres::SILENT;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            return {false, summary.FullReport()};
        }
        p[0] *= scale;
        residual = summary.final_cost;
        if (residuals.size() > 0) {
            std::vector<double> r_temp(G.size() + b.size());
            problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
            for (int i = 0; i < G.size(); i++)
                residuals[i] = r_temp[i];
            Eigen::ArrayXd as = (-model.sequence.TR / fixed[2]).exp();
            for (int i = 0; i < a.size(); i++) {
                residuals[i + G.size()] = as[i] - a[i];
            }
            for (int i = 0; i < b.size(); i++) {
                residuals[i + G.size() + a.size()] = r_temp[i + G.size()];
            }
        }
        iterations = summary.iterations.size();
        return {true, ""};
    }
};

int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates qMT parameters from ellipse parameters.\nInputs are G, "
                                "a, b.\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> G_path(parser, "G_FILE", "Input G file");
    args::Positional<std::string> a_path(parser, "a_FILE", "Input a file");
    args::Positional<std::string> b_path(parser, "b_FILE", "Input b file");

    QI_COMMON_ARGS;
    args::Flag debug(parser, "DEBUG", "Output debugging messages", {'d', "debug"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio)", {'b', "B1"});
    args::ValueFlag<std::string> f0(parser, "f0", "f0 map (in Hertz)", {'f', "f0"});
    args::ValueFlag<std::string> T2_f(parser, "T2f", "T2 Free map (for simulation only)", {"T2f"});
    args::ValueFlag<double>      T2_b_us(
        parser, "T2b", "T2 of bound pool (in microseconds, default 12)", {"T2b"}, 12);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::CheckPos(G_path);
    QI::CheckPos(a_path);
    QI::CheckPos(b_path);

    QI::Log(verbose, "Reading sequence information");
    json input = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);

    QI::SSFPMTSequence   ssfp(input.at("SSFPMT"));
    Eigen::ArrayXd const Wtot = pow(ssfp.FA / ssfp.pulse.p1, 2) * (ssfp.pulse.p2 / ssfp.Trf);
    QI_DB(ssfp.pulse.p1)
    QI_DB(ssfp.pulse.p2)
    QI_DBVEC(ssfp.FA);
    QI_DBVEC(ssfp.Trf);
    QI_DBVEC(Wtot);
    EMTModel model{ssfp, Wtot, T2_b_us.Get() * 1e-6};
    if (simulate) {
        QI::SimulateModel<EMTModel, true>(input,
                                          model,
                                          {f0.Get(), B1.Get(), T2_f.Get()},
                                          {G_path.Get(), a_path.Get(), b_path.Get()},
                                          verbose,
                                          simulate.Get());
    } else {
        // First calculate T2_f
        auto a_input = QI::ReadImage<QI::VectorVolumeF>(a_path.Get(), verbose);

        QI::VolumeF::Pointer T2_f_calc = QI::VolumeF::New();
        T2_f_calc->CopyInformation(a_input);
        T2_f_calc->SetRegions(a_input->GetBufferedRegion());
        T2_f_calc->Allocate(true);
        QI::Info(verbose, "Calculating T2_f");
        auto mt = itk::MultiThreaderBase::New();
        mt->SetNumberOfWorkUnits(threads.Get());
        mt->ParallelizeImageRegion<3>(
            a_input->GetBufferedRegion(),
            [&](const QI::VectorVolumeF::RegionType &region) {
                itk::ImageRegionConstIterator<QI::VectorVolumeF> a_it(a_input, region);
                itk::ImageRegionIterator<QI::VolumeF>            T2_f_it(T2_f_calc, region);
                for (a_it.GoToBegin(); !a_it.IsAtEnd(); ++a_it, ++T2_f_it) {
                    auto const as = Eigen::Map<const Eigen::ArrayXf>(a_it.Get().GetDataPointer(),
                                                                     a_it.Get().Size());
                    Eigen::ArrayXd T2_fs = (-ssfp.TR / as.cast<double>().log());
                    T2_f_it.Set(T2_fs.mean()); // Different TRs so have to average afterwards
                }
            },
            nullptr);

        EMTFit fit{model};
        auto   fit_filter = QI::ModelFitFilter<EMTFit>::New(&fit, verbose, resids, subregion.Get());
        fit_filter->ReadInputs(
            {G_path.Get(), a_path.Get(), b_path.Get()}, {f0.Get(), B1.Get(), ""}, mask.Get());
        fit_filter->SetFixed(2, T2_f_calc);
        fit_filter->Update();
        fit_filter->WriteOutputs(prefix.Get() + "EMT_");
        QI::WriteImage(T2_f_calc, prefix.Get() + "EMT_T2_f" + QI::OutExt(), verbose);
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
