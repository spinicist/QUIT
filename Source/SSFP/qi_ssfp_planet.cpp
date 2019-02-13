/*
 *  qdespot1.cpp - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <array>

#include <Eigen/Core>

#include "Args.h"
#include "ImageIO.h"
#include "ModelFitFilter.h"
#include "SSFPSequence.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

struct PLANETModel {
    using DataType      = double;
    using ParameterType = double;
    using SequenceType  = QI::SSFPSequence;

    static constexpr int NV = 3;
    static constexpr int ND = 0;
    static constexpr int NF = 1;

    static std::array<const std::string, NV> varying_names;
    static std::array<const std::string, NF> fixed_names;
    static const QI_ARRAYN(double, NF) fixed_defaults;
    const SequenceType &sequence;

    size_t num_outputs() const { return 3; }
    int    output_size(int /* Unused */) { return sequence.size(); }

    auto signals(const Eigen::ArrayBase<QI_ARRAYN(double, NV)> &v,
                 const QI_ARRAYN(double, NF) & f) const -> std::vector<QI_ARRAY(double)> {
        const double &                PD      = v[0];
        const double &                T1      = v[1];
        const double &                T2      = v[2];
        const double &                B1      = f[0];
        const double                  E1      = exp(-sequence.TR / T1);
        const double                  E2      = exp(-sequence.TR / T2);
        const Eigen::ArrayXd          alpha   = B1 * sequence.FA;
        const Eigen::ArrayXd          d       = (1. - E1 * E2 * E2 - (E1 - E2 * E2) * cos(alpha));
        const Eigen::ArrayXd          G       = PD * sqrt(E2) * (1 - E1) * sin(alpha) / d;
        const Eigen::ArrayXd          a       = Eigen::ArrayXd::Constant(sequence.size(), E2);
        const Eigen::ArrayXd          b       = E2 * (1. - E1) * (1. + cos(alpha)) / d;
        std::vector<QI_ARRAY(double)> outputs = {G, a, b};
        return outputs;
    }
};
std::array<const std::string, 3> PLANETModel::varying_names{{"PD"s, "T1"s, "T2"s}};
std::array<const std::string, 1> PLANETModel::fixed_names{{"B1"s}};
const QI_ARRAYN(double, 1) PLANETModel::fixed_defaults{1.0};

struct PLANETFit {
    static const bool Blocked = true;
    static const bool Indexed = false;
    using InputType           = double;
    using OutputType          = double;
    using ResidualType        = double;
    using FlagType            = int;
    using ModelType           = PLANETModel;
    ModelType model;

    int n_inputs() const { return 3; }
    int input_size(const int /* Unused */) const { return 1; }
    int n_fixed() const { return 1; }
    int n_outputs() const { return 3; }

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          const Eigen::ArrayXd &             fixed,
                          QI_ARRAYN(OutputType, PLANETModel::NV) & out,
                          ResidualType & /* Unused */,
                          std::vector<Eigen::ArrayXd> & /* Unused */,
                          FlagType & /* Unused */,
                          const int block) const {
        const double &G    = inputs[0][0];
        const double &a    = inputs[1][0];
        const double &b    = inputs[2][0];
        const double  b1   = fixed[0];
        const double  cosa = cos(b1 * model.sequence.FA(block));
        const double  sina = sin(b1 * model.sequence.FA(block));
        const double  T1   = -model.sequence.TR / log((a * (1. + cosa - a * b * cosa) - b) /
                                                   (a * (1. + cosa - a * b) - b * cosa));
        const double  T2   = -model.sequence.TR / log(a);
        const double  E1   = exp(-model.sequence.TR / T1);
        const double  E2   = a; // For simplicity copying formulas
        const double  PD =
            G * (1. - E1 * cosa - E2 * E2 * (E1 - cosa)) / (sqrt(E2) * (1. - E1) * sina);
        out[0] = PD;
        out[1] = T1;
        out[2] = T2;
        return std::make_tuple(true, "");
    }
};

int main(int argc, char **argv) {
    args::ArgumentParser parser(
        "Calculates T1&T2 from SSFP Ellipse Parameters.\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> G_path(parser, "G", "Ellipse parameter G");
    args::Positional<std::string> a_path(parser, "a", "Ellipse parameter a");
    args::Positional<std::string> b_path(parser, "b", "Ellipse parameter b");
    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> out_prefix(parser, "OUTPREFIX", "Add a prefix to output filenames",
                                            {'o', "out"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask",
                                      {'m', "mask"});
    args::ValueFlag<std::string> subregion(
        parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK",
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

    QI::Log(verbose, "Reading sequence information");
    rapidjson::Document input = seq_arg ? QI::ReadJSON(seq_arg.Get()) : QI::ReadJSON(std::cin);
    QI::SSFPSequence    ssfp(QI::GetMember(input, "SSFP"));
    PLANETModel         model{ssfp};
    if (simulate) {
        QI::SimulateModel<PLANETModel, true>(input, model, {B1.Get()},
                                             {G_path.Get(), a_path.Get(), b_path.Get()}, verbose,
                                             simulate.Get());
    } else {
        PLANETFit fit{model};
        auto fit_filter = QI::ModelFitFilter<PLANETFit>::New(&fit, verbose, false, subregion.Get());
        fit_filter->ReadInputs({G_path.Get(), a_path.Get(), b_path.Get()}, {B1.Get()}, mask.Get());
        fit_filter->SetBlocks(ssfp.size());
        fit_filter->Update();
        fit_filter->WriteOutputs(out_prefix.Get() + "PLANET_");
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
