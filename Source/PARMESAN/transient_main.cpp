/*
 *  transient_main.cpp
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
#include "FitScaledNumeric.h"
#include "ImageIO.h"
#include "Macro.h"
#include "ModelFitFilter.h"
#include "SimulateModel.h"
#include "Util.h"

#include "transient_b1_model.h"
#include "transient_model.h"
#include "transient_mt_model.h"

/*
 * Main
 */
int transient_main(args::Subparser &parser) {
    args::Positional<std::string> input_path(parser, "INPUT", "Input MUPA file");
    QI_COMMON_ARGS;
    args::Flag                   mt(parser, "MT", "Use MT model", {"mt"});
    args::ValueFlag<double>      T2_b(parser, "T2_b", "T2 of bound pool", {"T2b"}, 12e-6);
    args::ValueFlag<std::string> ls_arg(
        parser, "LINESHAPE", "Path to lineshape file", {"lineshape"});

    parser.Parse();

    QI::CheckPos(input_path);

    QI::Log(verbose, "Reading sequence parameters");
    json doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);

    RUFISSequence sequence(doc["MUPA"]);

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
            FitType fit{model};
            auto    fit_filter =
                QI::ModelFitFilter<FitType>::New(&fit, verbose, covar, resids, subregion.Get());
            fit_filter->ReadInputs({input_path.Get()}, fixed, mask.Get());
            fit_filter->Update();
            fit_filter->WriteOutputs(prefix.Get() + model_name);
        }
    };

    if (mt) {
        MUPAMTModel model{{}, sequence};
        process(model, "MUPAMT_", {});
    } else {
        MUPAB1Model model{{}, sequence};
        process(model, "MUPAB1_", {});
    }
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
