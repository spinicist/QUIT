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

int parmesan_fit(args::Subparser &parser) {
    args::Positional<std::string> json_path(parser, "JSON", "Parameter file");
    args::Positional<std::string> basis_path(parser, "BASIS", "Basis file");
    args::Positional<std::string> input_path(parser, "INPUT", "Input image file");
    QI_COMMON_ARGS;
    args::ValueFlag<std::string> mt(
        parser, "MT", "Use MT model (provide interp table for linearized GBM)", {"mt"});
    args::ValueFlag<double> T2_f(parser, "T2_f", "Fix T2_f", {"T2_f"}, 0.07);
    args::ValueFlag<double> T1_s(parser, "T1_s", "Fix T1_s", {"T1_s"}, 0.35);
    args::ValueFlag<double> T2_s(parser, "T2_s", "Fix T2_s", {"T2_s"}, 14e-6);
    args::ValueFlag<double> R_x(parser, "T2_f", "Fix R_x", {"R_x"}, 14);
    parser.Parse();
    QI::CheckPos(input_path);
    PrepZTESequence sequence(QI::ReadJSON(json_path.Get())["PrepZTE"]);

    auto process = [&](auto                                       model,
                       const std::string                         &model_name,
                       typename decltype(model)::FixedNames const fixed) {
        using FitType = QI::ScaledNumericDiffFit<decltype(model)>;
        FitType fit{model};
        auto    fit_filter =
            QI::ModelFitFilter<FitType>::New(&fit, covar, resids, threads.Get(), subregion.Get());
        fit_filter->ReadInputs({input_path.Get()}, fixed, mask.Get());
        fit_filter->SetRequestedRegionFromString(subregion.Get());
        fit_filter->Update();
        fmt::print(stderr, "Returned from update\n");
        fit_filter->WriteOutputs(prefix.Get() + model_name);
        fmt::print(stderr, "Wrote outputs\n");
    };

    if (mt) {
        json            lgbm = QI::ReadJSON(mt.Get());
        QI::RegularGrid R2sl(QI::ArrayFromJSON<double>(lgbm, "ω"),
                             QI::ArrayFromJSON<double>(lgbm, "τ"),
                             QI::MatrixFromJSON<double>(lgbm, "A"));
        if (T2_f || T1_s || T2_s || R_x) {
            QI::PrepQMTModel model{
                {}, sequence, R2sl, T2_f.Get(), T1_s.Get(), T2_s.Get(), R_x.Get()};
            process(model, "PARMESAN_", {});
        } else {
            QI::PrepQMTFullModel model{{}, sequence, R2sl};
            process(model, "PARMESAN_full_", {});
        }
    } else {
        PrepB1Model model{{}, sequence};
        process(model, "PARMESAN_", {});
    }
    QI::Info("Finished.");
    return EXIT_SUCCESS;
}

int parmesan_sim(args::Subparser &parser) {
    args::Positional<std::string>     json_path(parser, "JSON", "Parameter file");
    args::Positional<std::string>     basis_path(parser, "BASIS", "Basis file");
    args::Positional<std::string>     out_path(parser, "OUTPUT", "Simulation output file");
    args::PositionalList<std::string> varying_paths(parser, "INPUT", "Input parameter maps");
    QI_CORE_ARGS;
    args::ValueFlag<float> noise(parser, "NOISE", "Noise standard deviation", {'n', "noise"}, 0.f);
    args::ValueFlag<std::string> mt(parser, "MT", "MT model (interp table for Lin GBM)", {"mt"});
    args::ValueFlag<double>      T2_f(parser, "T2_f", "Fix T2_f", {"T2_f"}, 0.07);
    args::ValueFlag<double>      T1_s(parser, "T1_s", "Fix T1_s", {"T1_s"}, 0.35);
    args::ValueFlag<double>      T2_s(parser, "T2_s", "Fix T2_s", {"T2_s"}, 14e-6);
    args::ValueFlag<double>      R_x(parser, "T2_f", "Fix R_x", {"R_x"}, 14);
    parser.Parse();
    QI::Info("Reading sequence parameters");
    PrepZTESequence sequence(QI::ReadJSON(json_path.Get())["PrepZTE"]);

    auto process = [&](auto                                       model,
                       const std::string                         &model_name,
                       typename decltype(model)::FixedNames const fixed) {
        QI::SimulateModel2<decltype(model), false>(model,
                                                   varying_paths.Get(),
                                                   fixed,
                                                   {out_path.Get()},
                                                   mask.Get(),
                                                   noise.Get(),
                                                   threads.Get(),
                                                   subregion.Get());
    };

    if (mt) {
        json            lgbm = QI::ReadJSON(mt.Get());
        QI::RegularGrid R2sl(QI::ArrayFromJSON<double>(lgbm, "ω"),
                             QI::ArrayFromJSON<double>(lgbm, "τ"),
                             QI::MatrixFromJSON<double>(lgbm, "A"));
        if (T2_f || T1_s || T2_s || R_x) {
            QI::PrepQMTModel model{
                {}, sequence, R2sl, T2_f.Get(), T1_s.Get(), T2_s.Get(), R_x.Get()};
            process(model, "PARMESAN_", {});
        } else {
            QI::PrepQMTFullModel model{{}, sequence, R2sl};
            process(model, "PARMESAN_full_", {});
        }
    } else {
        PrepB1Model model{{}, sequence};
        process(model, "PARMESAN_", {});
    }
    QI::Info("Finished.");
    return EXIT_SUCCESS;
}
