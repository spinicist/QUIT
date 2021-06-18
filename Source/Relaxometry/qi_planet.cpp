/*
 *  qi_planet.cpp - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>
#include <array>

#include "Args.h"
#include "ImageIO.h"
#include "ModelFitFilter.h"
#include "SSFPSequence.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

struct PLANETModel : QI::Model<double, double, 3, 1, 3> {
    const QI::SSFPSequence &sequence;

    static std::array<const std::string, NV> varying_names;
    static std::array<const std::string, NF> fixed_names;
    static const QI_ARRAYN(double, NF) fixed_defaults;

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
    using RMSErrorType        = double;
    using FlagType            = int;
    using ModelType           = PLANETModel;
    ModelType model;

    int input_size(const int /* Unused */) const { return 1; }
    int n_outputs() const { return 3; }

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          const Eigen::ArrayXd &             fixed,
                          PLANETModel::VaryingArray &        out,
                          PLANETModel::CovarArray * /* Unused */,
                          RMSErrorType & /* Unused */,
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
        return {true, ""};
    }
};

int planet_main(args::Subparser &parser) {
    args::Positional<std::string> G_path(parser, "G", "Ellipse parameter G");
    args::Positional<std::string> a_path(parser, "a", "Ellipse parameter a");
    args::Positional<std::string> b_path(parser, "b", "Ellipse parameter b");
    args::ValueFlag<std::string>  B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    QI_COMMON_ARGS;
    parser.Parse();
    QI::CheckPos(G_path);
    QI::CheckPos(a_path);
    QI::CheckPos(b_path);

    QI::Log(verbose, "Reading sequence information");
    json             input = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    QI::SSFPSequence ssfp(input.at("SSFP"));
    PLANETModel      model{{}, ssfp};
    if (simulate) {
        QI::SimulateModel<PLANETModel, true>(input,
                                             model,
                                             {B1.Get()},
                                             {G_path.Get(), a_path.Get(), b_path.Get()},
                                             mask.Get(),
                                             verbose,
                                             simulate.Get(),
                                             subregion.Get());
    } else {
        PLANETFit fit{model};
        auto      fit_filter =
            QI::ModelFitFilter<PLANETFit>::New(
                &fit, verbose, false, false, threads.Get(), subregion.Get());
        fit_filter->ReadInputs({G_path.Get(), a_path.Get(), b_path.Get()}, {B1.Get()}, mask.Get());
        fit_filter->SetBlocks(ssfp.size());
        fit_filter->Update();
        fit_filter->WriteOutputs(prefix.Get() + "PLANET_");
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
