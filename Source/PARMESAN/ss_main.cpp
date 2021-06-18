/*
 *  ss_main.cpp
 *
 *  Copyright (c) 2020 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>

// #define QI_DEBUG_BUILD 1

#include "Args.h"
#include "FitScaledNumeric.h"
#include "ImageIO.h"
#include "Macro.h"
#include "ModelFitFilter.h"
#include "SimulateModel.h"
#include "Util.h"

#include "ss_T2.h"
#include "ss_model.h"
#include "ss_mt.h"

/*
 * Main
 */
int ss_main(args::Subparser &parser) {
    args::Positional<std::string> input_path(parser, "INPUT", "Input MUPA file");
    QI_COMMON_ARGS;
    args::Flag                   T2(parser, "T2", "Fit T2 model", {"T2"});
    args::Flag                   MT(parser, "MT", "Fit MT model", {"MT"});
    args::ValueFlag<std::string> ls_arg(
        parser, "LINESHAPE", "Path to lineshape file", {"lineshape"});
    parser.Parse();
    QI::CheckPos(input_path);
    QI::Log(verbose, "Reading sequence parameters");
    json       doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    SSSequence sequence(doc);

    auto process = [&](auto                                       model,
                       const std::string &                        model_name,
                       typename decltype(model)::FixedNames const fixed) {
        if (simulate) {
            QI::SimulateModel<decltype(model), false>(doc,
                                                      model,
                                                      fixed,
                                                      {input_path.Get()},
                                                      mask.Get(),
                                                      verbose,
                                                      simulate.Get(),
                                                      subregion.Get());
        } else {
            using FitType = QI::ScaledNumericDiffFit<decltype(model), decltype(model)::NS>;
            FitType fit(model);

            auto fit_filter =
                QI::ModelFitFilter<FitType>::New(
                    &fit, verbose, covar, resids, threads.Get(), subregion.Get());
            fit_filter->ReadInputs({input_path.Get()}, fixed, mask.Get());
            fit_filter->Update();
            fit_filter->WriteOutputs(prefix.Get() + model_name);
        }
    };

    if (MT) {
        QI::Log(verbose, "Using MT model");
        QI::Log(verbose, "Reading lineshape file: {}", ls_arg.Get());
        json        ls_file = QI::ReadJSON(ls_arg.Get());
        SS_MT_Model model{{}, sequence, ls_file.at("lineshape").get<QI::InterpLineshape>()};
        process(model, "SS_", {});
    } else if (T2) {
        QI::Log(verbose, "Using T2 model");
        SS_T1T2_Model model{{}, sequence};
        process(model, "SS_", {});
    } else {
        QI::Log(verbose, "Using T1 model");
        SS_T1_Model model{{}, sequence}; //, ls_file.at("lineshape").get<QI::InterpLineshape>()};
        process(model, "SS_", {});
    }

    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
