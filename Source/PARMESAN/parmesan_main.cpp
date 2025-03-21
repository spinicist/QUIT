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
#include "FitScaledAuto.h"
#include "FitScaledNumeric.h"
#include "ImageIO.h"
#include "Macro.h"
#include "ModelFitFilter.h"
#include "SimulateModel.h"
#include "Util.h"

#include "prep_b1_model.h"
#include "prep_qmt_model.h"

/*
 * Main
 */
int parmesan_main(args::Subparser &parser) {
    args::Positional<std::string> input_path(parser, "INPUT", "Input file");
    args::Positional<std::string> basis_path(parser, "BASIS", "Basis file");
    QI_COMMON_ARGS;
    args::ValueFlag<std::string> mt(parser, "MT", "Use MT model (provide interp table for linearized GBM)", {"mt"});
    parser.Parse();
    QI::CheckPos(input_path);
    QI::Log(verbose, "Reading sequence parameters");
    json doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);

    PrepZTESequence sequence(doc["PrepZTE"]);

    auto process = [&](auto                                       model,
                       const std::string                         &model_name,
                       typename decltype(model)::FixedNames const fixed) {
        if (simulate) {
            QI::SimulateModel<decltype(model), false>(doc,
                                                      model,
                                                      fixed,
                                                      {input_path.Get()},
                                                      mask.Get(),
                                                      verbose,
                                                      simulate.Get(),
                                                      threads.Get(),
                                                      subregion.Get());
        } else {
            using FitType = QI::ScaledNumericDiffFit<decltype(model)>;
            FitType fit{model};
            auto    fit_filter = QI::ModelFitFilter<FitType>::New(
                &fit, verbose, covar, resids, threads.Get(), subregion.Get());
            fit_filter->ReadInputs({input_path.Get()}, fixed, mask.Get());
            fit_filter->Update();
            fit_filter->WriteOutputs(prefix.Get() + model_name);
        }
    };

    if (mt) {
        json lgbm = QI::ReadJSON(mt.Get());
        QI::RegularGrid R2sl(QI::ArrayFromJSON<double>(lgbm, "tau"), QI::ArrayFromJSON<double>(lgbm, "omega"), QI::MatrixFromJSON<double>(lgbm, "R2sl"));
        QI::PrepQMTModel model{{}, sequence, R2sl};
        process(model, "PARMESAN_", {});
    } else {
        PrepB1Model model{{}, sequence};
        process(model, "PARMESAN_", {});
    }
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
