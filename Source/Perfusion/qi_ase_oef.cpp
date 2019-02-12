/*
 *  qi_ase_oef.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Dense>

// #define QI_DEBUG_BUILD

#include "Args.h"
#include "FitFunction.h"
#include "ImageIO.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "MultiEchoSequence.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

constexpr double kappa      = 0.03;                // Conversion factor
constexpr double gyro_gamma = 2 * M_PI * 42.577e6; // Gyromagnetic Ratio
constexpr double delta_X0   = 0.264e-6;            // Susc diff oxy and fully de-oxy blood
constexpr double Hct        = 0.34;
constexpr double Hb         = Hct / kappa;

struct ASEModel {
    using SequenceType  = QI::MultiEchoFlexSequence;
    using DataType      = double;
    using ParameterType = double;

    const SequenceType &sequence;
    const double        B0;

    static constexpr int NV = 4;
    static constexpr int ND = 3;
    static constexpr int NF = 1;
    using VaryingArray      = QI_ARRAYN(ParameterType, NV);
    using DerivedArray      = QI_ARRAYN(ParameterType, ND);
    using FixedArray        = QI_ARRAYN(ParameterType, NF);
    const std::array<const std::string, NV> varying_names{{"S0"s, "dT"s, "R2p"s, "DBV"s}};
    const std::array<const std::string, ND> derived_names{{"Tc"s, "OEF"s, "dHb"s}};
    const std::array<const std::string, NF> fixed_names{{}};
    const FixedArray                        fixed_defaults{};

    const VaryingArray start{0.98, 0., 1.0, 0.015};
    const VaryingArray bounds_lo{0.1, -0.1, 0.25, 0.001};
    const VaryingArray bounds_hi{2., 0.1, 10.0, 0.05};

    template <typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &varying, const FixedArray & /* Unused */) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T      = typename Derived::Scalar;
        const T &S0  = varying[0];
        const T &dT  = varying[1];
        const T &R2p = varying[2];
        const T &DBV = varying[3];

        const auto dw = R2p / DBV;
        const auto Tc = 1.5 * (1. / dw);

        const auto aTE = (sequence.TE + dT).abs();
        QI_ARRAY(T) S(sequence.size());
        for (int i = 0; i < sequence.size(); i++) {
            const auto tau = aTE(i);
            if (tau < Tc) {
                S(i) = S0 * exp(-(0.3) * (R2p * tau) * (R2p * tau) / DBV);
            } else {
                S(i) = S0 * exp(-tau * R2p + DBV);
            }
        }
        return S;
    }

    void derived(const VaryingArray &varying, const FixedArray & /* Unused */,
                 DerivedArray &      derived) const {
        const auto &R2p = varying[2];
        const auto &DBV = varying[3];

        const auto dw  = R2p / DBV;
        const auto Tc  = 1.5 * (1. / dw);
        const auto OEF = dw / ((4. * M_PI / 3.) * gyro_gamma * B0 * delta_X0 * Hct);
        const auto dHb = OEF * Hb;

        derived[0] = Tc;
        derived[1] = OEF;
        derived[2] = dHb;
    }
};
using ASEFit = QI::ScaledNLLSFitFunction<ASEModel>;

struct ASEFixDBVModel {
    using SequenceType  = QI::MultiEchoFlexSequence;
    using DataType      = double;
    using ParameterType = double;

    static constexpr int NV = 3;
    static constexpr int ND = 3;
    static constexpr int NF = 1;
    using VaryingArray      = QI_ARRAYN(ParameterType, NV);
    using DerivedArray      = QI_ARRAYN(ParameterType, ND);
    using FixedArray        = QI_ARRAYN(ParameterType, NF);

    const SequenceType &sequence;
    const double        B0, DBV;

    const std::array<const std::string, NV> varying_names{{"S0"s, "dT"s, "R2p"s}};
    const std::array<const std::string, ND> derived_names{{"Tc"s, "OEF"s, "dHb"s}};
    const std::array<const std::string, NF> fixed_names{{}};
    const FixedArray                        fixed_defaults{};

    const VaryingArray start{0.98, 0., 1.0};
    const VaryingArray bounds_lo{0.1, -0.1, 0.25};
    const VaryingArray bounds_hi{2., 0.1, 10.0};

    template <typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &varying, const FixedArray & /* Unused */) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T        = typename Derived::Scalar;
        const T &  S0  = varying[0];
        const T &  dT  = varying[1];
        const T &  R2p = varying[2];
        const auto dw  = R2p / DBV;
        const auto Tc  = 1.5 * (1. / dw);

        const auto aTE = (sequence.TE + dT).abs();
        QI_ARRAY(T) S(sequence.size());
        for (int i = 0; i < sequence.size(); i++) {
            const auto tau = aTE(i);
            if (tau < Tc) {
                S(i) = S0 * exp(-(0.3) * (R2p * tau) * (R2p * tau) / DBV);
            } else {
                S(i) = S0 * exp(-tau * R2p + DBV);
            }
        }
        return S;
    }

    void derived(const VaryingArray &varying, const FixedArray & /* Unused */,
                 DerivedArray &      derived) const {
        const auto &R2p = varying[2];
        const auto  dw  = R2p / DBV;
        const auto  Tc  = 1.5 * (1. / dw);
        const auto  OEF = dw / ((4. * M_PI / 3.) * gyro_gamma * B0 * delta_X0 * Hct);
        const auto  dHb = OEF * Hb;

        derived[0] = Tc;
        derived[1] = OEF;
        derived[2] = dHb;
    }
};
using ASEFixDBVFit = QI::ScaledNLLSFitFunction<ASEFixDBVModel>;

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser(
        "Calculates the OEF from ASE data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "ASE_FILE", "Input ASE file");
    args::HelpFlag                help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filename",
                                        {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask",
                                      {'m', "mask"});
    args::ValueFlag<double> B0(parser, "B0", "Field-strength (Tesla), default 3", {'B', "B0"}, 3.0);
    args::ValueFlag<double> DBV(parser, "DBV", "Fix DBV and only fit R2'", {'d', "DBV"}, 0.0);
    args::ValueFlag<std::string> gradz(
        parser, "GRADZ", "Gradient of field-map in z-direction for MFG correction", {'z', "gradz"});
    args::ValueFlag<double> slice_arg(
        parser, "SLICE THICKNESS",
        "Slice-thickness for MFG calculation (useful if there was a slice gap)", {'s', "slice"});
    args::ValueFlag<std::string> subregion(
        parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK",
        {'s', "subregion"});
    args::ValueFlag<std::string> json_file(parser, "FILE",
                                           "Read JSON input from file instead of stdin", {"file"});
    args::ValueFlag<float>       simulate(
        parser, "SIMULATE", "Simulate sequence instead of fitting model (argument is noise level)",
        {"simulate"}, 0.0);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    rapidjson::Document json = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    QI::MultiEchoFlexSequence sequence(json["MultiEchoFlex"]);

    if (simulate) {
        ASEModel model{sequence, B0.Get()};
        QI::SimulateModel<ASEModel, false>(json, model, {gradz.Get()}, {input_path.Get()}, verbose,
                                           simulate.Get());
    } else {
        if (DBV) {
            ASEFixDBVModel model{sequence, B0.Get(), DBV.Get()};
            ASEFixDBVFit   fit(model);
            auto           fit_filter = QI::ModelFitFilter<ASEFixDBVFit>::New(&fit, verbose, false);
            QI::Log(verbose, "Reading ASE data from: {}", QI::CheckPos(input_path));
            auto input = QI::ReadVectorImage(QI::CheckPos(input_path));
            // QI::VolumeF::SpacingType  vox_size = input->GetSpacing();
            // if (slice_arg) {
            //     vox_size[2] = slice_arg.Get();
            // }
            fit_filter->SetInput(0, input);
            if (mask)
                fit_filter->SetMask(QI::ReadImage(mask.Get(), verbose));
            if (gradz)
                fit_filter->SetFixed(2, QI::ReadImage(gradz.Get(), verbose));
            if (subregion) {
                fit_filter->SetSubregion(QI::RegionArg(args::get(subregion)));
            }
            fit_filter->Update();

            const std::string outPrefix = outarg ? outarg.Get() : QI::Basename(input_path.Get());
            for (size_t i = 0; i < model.NV; i++) {
                const std::string fname = outPrefix + "_" + model.varying_names[i] + QI::OutExt();
                QI::WriteImage(fit_filter->GetOutput(i), fname, verbose);
            }
            for (size_t i = 0; i < model.ND; i++) {
                const std::string fname = outPrefix + "_" + model.derived_names[i] + QI::OutExt();
                QI::WriteImage(fit_filter->GetDerivedOutput(i), fname, verbose);
            }
        } else {
            ASEModel model{sequence, B0.Get()};
            ASEFit   fit(model);
            auto     fit_filter = QI::ModelFitFilter<ASEFit>::New(&fit, verbose, false);
            QI::Log(verbose, "Reading ASE data from: {}", QI::CheckPos(input_path));
            auto input = QI::ReadVectorImage(QI::CheckPos(input_path));
            // QI::VolumeF::SpacingType  vox_size = input->GetSpacing();
            // if (slice_arg) {
            //     vox_size[2] = slice_arg.Get();
            // }
            fit_filter->SetInput(0, input);
            if (mask)
                fit_filter->SetMask(QI::ReadImage(mask.Get(), verbose));
            if (gradz)
                fit_filter->SetFixed(2, QI::ReadImage(gradz.Get(), verbose));
            if (subregion) {
                fit_filter->SetSubregion(QI::RegionArg(args::get(subregion)));
            }
            fit_filter->Update();

            const std::string outPrefix = outarg ? outarg.Get() : QI::Basename(input_path.Get());
            for (size_t i = 0; i < model.NV; i++) {
                const std::string fname = outPrefix + "_" + model.varying_names[i] + QI::OutExt();
                QI::WriteImage(fit_filter->GetOutput(i), fname, verbose);
            }
            for (size_t i = 0; i < model.ND; i++) {
                const std::string fname = outPrefix + "_" + model.derived_names[i] + QI::OutExt();
                QI::WriteImage(fit_filter->GetDerivedOutput(i), fname, verbose);
            }
            return EXIT_SUCCESS;
        }
    }
}
