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
#include "prep_qmt_Rx.h"
#include "prep_qmt_k.h"

int parmesan_fit(args::Subparser &parser) {
    args::Positional<std::string> json_path(parser, "JSON", "Parameter file");
    args::Positional<std::string> input_path(parser, "INPUT", "Input image file");
    QI_COMMON_ARGS;
    Parse(parser);
    QI::CheckPos(input_path);
    PrepZTESequence sequence(QI::ReadJSON(json_path.Get())["PrepZTE"]);

    PrepB1Model model{{}, sequence};
    using FitType = QI::ScaledNumericDiffFit<decltype(model)>;
    FitType fit{model};
    auto    fit_filter =
        QI::ModelFitFilter<FitType>::New(&fit, covar, resids, threads.Get(), subregion.Get());
    fit_filter->ReadInputs({input_path.Get()}, {}, mask.Get());
    fit_filter->SetRequestedRegionFromString(subregion.Get());
    fit_filter->Update();
    fit_filter->WriteOutputs(prefix.Get() + "PARMESAN_");

    QI::Info("Finished.");
    return EXIT_SUCCESS;
}

int parmesan_sim(args::Subparser &parser) {
    args::Positional<std::string>     json_path(parser, "JSON", "Parameter file");
    args::Positional<std::string>     out_path(parser, "OUTPUT", "Simulation output file");
    args::PositionalList<std::string> varying_paths(parser, "INPUT", "Input parameter maps");
    QI_CORE_ARGS;
    args::ValueFlag<float> noise(parser, "NOISE", "Noise standard deviation", {'n', "noise"}, 0.f);
    parser.Parse();
    QI::Info("Reading sequence parameters");
    PrepZTESequence sequence(QI::ReadJSON(json_path.Get())["PrepZTE"]);
    PrepB1Model     model{{}, sequence};
    QI::SimulateModel2<decltype(model), false>(model,
                                               varying_paths.Get(),
                                               {},
                                               {out_path.Get()},
                                               mask.Get(),
                                               noise.Get(),
                                               threads.Get(),
                                               subregion.Get());
    QI::Info("Finished.");
    return EXIT_SUCCESS;
}

int parmesan_qmt_fit(args::Subparser &parser) {
    args::Positional<std::string> json_path(parser, "JSON", "Parameter file");
    args::Positional<std::string> input_path(parser, "INPUT", "Input image file");
    args::Positional<std::string> sl(parser, "MT", "Interp table for linearized GBM", {"sl"});
    QI_COMMON_ARGS;

    args::ValueFlag<double> R_x(parser, "T2_f", "Fix R_x", {"R_x"}, 14);
    args::ValueFlag<double> k(parser, "T2_f", "Fix k", {"k"}, 1.4);
    args::ValueFlag<double> T2_f(parser, "T2_f", "Fix T2_f", {"T2_f"}, 0.1);
    args::ValueFlag<double> T1_s(parser, "T1_s", "Fix T1_s", {"T1_s"}, 0.35);
    args::ValueFlag<double> T2_s(parser, "T2_s", "Fix T2_s", {"T2_s"}, 14e-6);
    Parse(parser);
    QI::CheckPos(input_path);
    PrepZTESequence sequence(QI::ReadJSON(json_path.Get())["PrepZTE"]);
    json            lgbm = QI::ReadJSON(sl.Get());
    QI::RegularGrid R2sl(QI::ArrayFromJSON<double>(lgbm, "ω"),
                         QI::ArrayFromJSON<double>(lgbm, "τ"),
                         QI::MatrixFromJSON<double>(lgbm, "A"));
    auto            process = [&](auto                                       model,
                       const std::string                         &model_name,
                       typename decltype(model)::FixedNames const fixed) {
        using FitType = QI::ScaledNumericDiffFit<decltype(model)>;
        FitType fit{model};
        auto    fit_filter =
            QI::ModelFitFilter<FitType>::New(&fit, covar, resids, threads.Get(), subregion.Get());
        fit_filter->ReadInputs({input_path.Get()}, fixed, mask.Get());
        fit_filter->SetRequestedRegionFromString(subregion.Get());
        fit_filter->Update();
        fit_filter->WriteOutputs(prefix.Get() + model_name);
    };

    if (R_x) {
        QI::PrepQMTRx model{{}, sequence, R2sl, T2_f.Get(), T1_s.Get(), T2_s.Get(), R_x.Get()};
        process(model, "PARMESAN_R_x_", {});
    } else if (k) {
        QI::PrepQMTk model{{}, sequence, R2sl, T2_f.Get(), T1_s.Get(), T2_s.Get(), k.Get()};
        process(model, "PARMESAN_k_", {});
    } else {
        QI::PrepQMTRxFull model{{}, sequence, R2sl};
        process(model, "PARMESAN_qMT_", {});
    }
    QI::Info("Finished.");
    return EXIT_SUCCESS;
}

int parmesan_qmt_sim(args::Subparser &parser) {
    args::Positional<std::string> json_path(parser, "JSON", "Parameter file");
    args::Positional<std::string> out_path(parser, "OUTPUT", "Simulation output file");
    args::Positional<std::string> mt(parser, "MT", "MT model (interp table for Lin GBM)", {"mt"});
    args::PositionalList<std::string> varying_paths(parser, "INPUT", "Input parameter maps");

    QI_CORE_ARGS;
    args::ValueFlag<float> noise(parser, "NOISE", "Noise standard deviation", {'n', "noise"}, 0.f);

    args::ValueFlag<std::string> T2_f(parser, "T2_f", "T2_f map", {"T2_f"});
    args::ValueFlag<std::string> T1_s(parser, "T1_s", "T1_s map", {"T1_s"});
    args::ValueFlag<std::string> T2_s(parser, "T2_s", "T2_s map", {"T2_s"});
    args::ValueFlag<std::string> R_x(parser, "R_x", "R_x map", {"R_x"});
    args::ValueFlag<std::string> k(parser, "k", "k map", {"k"});
    parser.Parse();
    QI::Info("Reading sequence parameters");
    PrepZTESequence sequence(QI::ReadJSON(json_path.Get())["PrepZTE"]);
    json            lgbm = QI::ReadJSON(mt.Get());
    QI::RegularGrid R2sl(QI::ArrayFromJSON<double>(lgbm, "ω"),
                         QI::ArrayFromJSON<double>(lgbm, "τ"),
                         QI::MatrixFromJSON<double>(lgbm, "A"));

    auto process = [&](auto model, typename decltype(model)::FixedNames const fixed) {
        QI::SimulateModel2<decltype(model), false>(model,
                                                   varying_paths.Get(),
                                                   fixed,
                                                   {out_path.Get()},
                                                   mask.Get(),
                                                   noise.Get(),
                                                   threads.Get(),
                                                   subregion.Get());
    };

    if (R_x) {
        QI::PrepQMTRx model{{}, sequence, R2sl};
        process(model, {T2_f.Get(), T1_s.Get(), T2_s.Get(), R_x.Get()});
    } else if (k) {
        QI::PrepQMTk model{{}, sequence, R2sl};
        process(model, {T2_f.Get(), T1_s.Get(), T2_s.Get(), k.Get()});
    } else {
        QI::PrepQMTRxFull model{{}, sequence, R2sl};
        process(model, {});
    }
    QI::Info("Finished.");
    return EXIT_SUCCESS;
}
