/*
 *  mupa_main.cpp
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>

// #define QI_DEBUG_BUILD 1

#include "Args.h"
#include "ImageIO.h"
#include "Macro.h"
#include "ModelFitFilter.h"
#include "SimulateModel.h"
#include "Util.h"

#include "mupa_model_1c.h"
#include "mupa_model_b1.h"
#include "mupa_model_mt.h"
#include "mupa_model_mtb1.h"

/*
 * Main
 */
int mupa_main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates parametric maps from MUPA data "
                                "data.\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> input_path(parser, "INPUT", "Input MUPA file");

    QI_COMMON_ARGS;

    args::ValueFlag<std::string> B1(parser, "B1", "Fix B1", {"B1"});
    args::Flag                   mt(parser, "MT", "Use MT model", {"mt"});
    args::ValueFlag<double>      T2_b(parser, "T2_b", "T2 of bound pool", {"T2b"}, 12e-6);
    args::ValueFlag<std::string> ls_arg(
        parser, "LINESHAPE", "Path to lineshape file", {"lineshape"});

    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI::CheckPos(input_path);

    QI::Log(verbose, "Reading sequence parameters");
    json doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);

    MUPASequence sequence(doc["MUPA"]);
    auto         process = [&](auto                                       model,
                       const std::string &                        model_name,
                       typename decltype(model)::FixedNames const fixed) {
        if (simulate) {
            QI::SimulateModel<decltype(model), false>(
                doc, model, fixed, {input_path.Get()}, verbose, simulate.Get());
        } else {
            using FitType = QI::ScaledNumericDiffFit<decltype(model), 2>;
            FitType fit{model};
            auto    fit_filter =
                QI::ModelFitFilter<FitType>::New(&fit, verbose, covar, resids, subregion.Get());
            fit_filter->ReadInputs({input_path.Get()}, fixed, mask.Get());
            fit_filter->Update();
            fit_filter->WriteOutputs(prefix.Get() + model_name);
        }
    };

    if (mt && B1) {
        QI::Log(verbose, "Reading lineshape file: {}", ls_arg.Get());
        json          ls_file = QI::ReadJSON(ls_arg.Get());
        MUPAMTB1Model model{{}, sequence, ls_file.at("lineshape").get<QI::InterpLineshape>()};
        process(model, "MUPAMTB1_", {B1.Get()});
    } else if (mt) {
        QI::Log(verbose, "Reading lineshape file: {}", ls_arg.Get());
        json        ls_file = QI::ReadJSON(ls_arg.Get());
        MUPAMTModel model{{}, sequence, ls_file.at("lineshape").get<QI::InterpLineshape>()};
        process(model, "MUPAMT_", {});
    } else {
        MUPAB1Model model{{}, sequence};
        process(model, "MUPAB1_", {});
    }
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
