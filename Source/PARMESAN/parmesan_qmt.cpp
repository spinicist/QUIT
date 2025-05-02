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

// #define QI_DEBUG_BUILD 1

#include "Args.h"
#include "FitScaledAuto.h"
#include "FitScaledNumeric.h"
#include "ImageIO.h"
#include "Macro.h"
#include "ModelFitFilter.h"
#include "SimulateModel.h"
#include "Util.h"

#include "Basis.h"
#include "ParameterGrid.h"
#include "prep_qmt.h"

#include <Eigen/Core>
#include <Eigen/SVD>
#include <functional>

int parmesan_qmt_fit(args::Subparser &parser) {
    args::Positional<std::string> json_path(parser, "JSON", "Parameter file");
    args::Positional<std::string> input_path(parser, "INPUT", "Input image file");
    args::Positional<std::string> sl(parser, "SL", "Super-Lorentzian Interp Table", {"sl"});
    QI_COMMON_ARGS;

    args::ValueFlag<std::string> bpath(parser, "BASIS", "Path to basis JSON file", {'b', "basis"});
    args::ValueFlag<double>      R_x(parser, "T2_f", "Fix R_x", {"R_x"}, 14);
    args::ValueFlag<double>      k(parser, "T2_f", "Fix k", {"k"}, 1.4);
    args::ValueFlag<double>      T2_f(parser, "T2_f", "Fix T2_f", {"T2_f"}, 0.1);
    args::ValueFlag<double>      T1_s(parser, "T1_s", "Fix T1_s", {"T1_s"}, 0.35);
    args::ValueFlag<double>      T2_s(parser, "T2_s", "Fix T2_s", {"T2_s"}, 14e-6);
    args::Flag                   kfull(parser, "k", "Use full k model", {"kfull"});
    Parse(parser);
    QI::CheckPos(input_path);
    PrepSequence   sequence(QI::ReadJSON(json_path.Get())["PrepZTE"]);
    json           lgbm = QI::ReadJSON(sl.Get());
    QI::InterpGrid R2sl(QI::ArrayFromJSON<double>(lgbm, "ω"),
                        QI::ArrayFromJSON<double>(lgbm, "τ"),
                        QI::MatrixFromJSON<double>(lgbm, "A"));
    auto           process = [&](auto                                       model,
                       const std::string                         &model_name,
                       typename decltype(model)::FixedNames const fixed) {
        using FitType = QI::ScaledNumericDiffFit<decltype(model)>;
        FitType fit{model};
        auto    fit_filter = QI::ModelFitFilter<FitType>::New(&fit, covar, resids, subregion.Get());
        fit_filter->ReadInputs({input_path.Get()}, fixed, mask.Get());
        fit_filter->SetRequestedRegionFromString(subregion.Get());
        fit_filter->Update();
        fit_filter->WriteOutputs(prefix.Get() + model_name);
    };

    if (R_x) {
        QI::PrepQMTRx model{{},
                            sequence,
                            R2sl,
                            ReadBasis(bpath.Get()),
                            T2_f.Get(),
                            T1_s.Get(),
                            T2_s.Get(),
                            R_x.Get()};
        process(model, "PARMESAN_R_x_", {});
    } else if (k) {
        QI::PrepQMTk model{{},
                           sequence,
                           R2sl,
                           ReadBasis(bpath.Get()),
                           T2_f.Get(),
                           T1_s.Get(),
                           T2_s.Get(),
                           k.Get()};
        process(model, "PARMESAN_k_", {});
    } else if (kfull) {
        QI::PrepQMTkFull model{{}, sequence, R2sl};
        process(model, "PARMESAN_qMT_", {});
    } else {
        QI::PrepQMTRxFull model{{}, sequence, R2sl};
        process(model, "PARMESAN_qMT_", {});
    }
    QI::Info("Finished.");
    return EXIT_SUCCESS;
}

int parmesan_qmt_sim(args::Subparser &parser) {
    args::Positional<std::string>     json_path(parser, "JSON", "Parameter file");
    args::Positional<std::string>     out_path(parser, "OUTPUT", "Simulation output file");
    args::Positional<std::string>     sl(parser, "SL", "Super-Lorentzian Interp Table", {"sl"});
    args::PositionalList<std::string> varying_paths(parser, "INPUT", "Input parameter maps");

    QI_CORE_ARGS;
    args::ValueFlag<float> noise(parser, "NOISE", "Noise standard deviation", {'n', "noise"}, 0.f);
    args::ValueFlag<std::string> bpath(parser, "BASIS", "Path to basis JSON file", {'b', "basis"});
    args::ValueFlag<std::string> T2_f(parser, "T2_f", "T2_f map", {"T2_f"});
    args::ValueFlag<std::string> T1_s(parser, "T1_s", "T1_s map", {"T1_s"});
    args::ValueFlag<std::string> T2_s(parser, "T2_s", "T2_s map", {"T2_s"});
    args::ValueFlag<std::string> R_x(parser, "R_x", "R_x map", {"R_x"});
    args::ValueFlag<std::string> k(parser, "k", "k map", {"k"});
    args::Flag                   kfull(parser, "k", "Use full k model", {"kfull"});
    Parse(parser);
    QI::Info("Reading sequence parameters");
    PrepSequence   sequence(QI::ReadJSON(json_path.Get())["PrepZTE"]);
    json           lgbm = QI::ReadJSON(sl.Get());
    QI::InterpGrid R2sl(QI::ArrayFromJSON<double>(lgbm, "ω"),
                        QI::ArrayFromJSON<double>(lgbm, "τ"),
                        QI::MatrixFromJSON<double>(lgbm, "A"));
    auto const     basis = ReadBasis(bpath.Get());

    auto process = [&](auto model, typename decltype(model)::FixedNames const fixed) {
        QI::SimulateModel2<decltype(model), false>(model,
                                                   varying_paths.Get(),
                                                   fixed,
                                                   {out_path.Get()},
                                                   mask.Get(),
                                                   noise.Get(),
                                                   subregion.Get());
    };

    if (R_x) {
        QI::PrepQMTRx model{{}, sequence, R2sl, basis};
        process(model, {T2_f.Get(), T1_s.Get(), T2_s.Get(), R_x.Get()});
    } else if (k) {
        QI::PrepQMTk model{{}, sequence, R2sl, basis};
        process(model, {T2_f.Get(), T1_s.Get(), T2_s.Get(), k.Get()});
    } else if (kfull) {
        QI::PrepQMTkFull model{{}, sequence, R2sl};
        process(model, {});
    } else {
        QI::PrepQMTRxFull model{{}, sequence, R2sl};
        process(model, {});
    }
    QI::Info("Finished.");
    return EXIT_SUCCESS;
}

int parmesan_qmt_basis(args::Subparser &parser) {
    args::Positional<std::string> json_path(parser, "JSON", "Parameter JSON file");
    args::Positional<std::string> out_path(parser, "OUTPUT", "Basis JSON file");
    args::Positional<std::string> sl(parser, "SL", "Super-Lorentzian Interp Table", {"sl"});
    args::ValueFlag<int>          B(parser, "B", "Basis size (8)", {'b', "B"}, 8);
    args::ValueFlag<int>          N(parser, "N", "Use N random samples", {'n', "N"}, 16384);
    args::Flag                    k(parser, "k", "Use k instead of R_x", {'k', "k"});

    Parse(parser);
    QI::CheckPos(json_path);
    QI::CheckPos(out_path);
    QI::Info("Reading sequence parameters");
    PrepSequence    sequence(QI::ReadJSON(json_path.Get())["PrepZTE"]);
    json            lgbm = QI::ReadJSON(sl.Get());
    QI::InterpGrid  R2sl(QI::ArrayFromJSON<double>(lgbm, "ω"),
                        QI::ArrayFromJSON<double>(lgbm, "τ"),
                        QI::MatrixFromJSON<double>(lgbm, "A"));
    Eigen::MatrixXd signals;
    Eigen::VectorXd ref;
    if (k) {
        QI::PrepQMTkFull model{{}, sequence, R2sl};
        auto const       pars = N ? QI::RandomPars(model.lo, model.hi, N.Get()) :
                                    QI::RegularPars(model.lo,
                                              model.hi,
                                              Eigen::Array<int, 9, 1>{1, 5, 5, 5, 5, 5, 5, 5, 5});
        signals.resize(pars.cols(), sequence.size());
        for (int ii = 0; ii < pars.cols(); ii++) {
            signals.row(ii) = model.signal(pars.col(ii), Eigen::ArrayXd{});
        }
        ref = signals.row(pars.cols() / 2);
    } else {
        QI::PrepQMTRxFull model{{}, sequence, R2sl};
        auto const        pars = N ? QI::RandomPars(model.lo, model.hi, N.Get()) :
                                     QI::RegularPars(model.lo,
                                              model.hi,
                                              Eigen::Array<int, 9, 1>{1, 5, 5, 5, 5, 5, 5, 5, 5});
        signals.resize(pars.cols(), sequence.size());
        for (int ii = 0; ii < pars.cols(); ii++) {
            auto sig = model.signal(pars.col(ii), Eigen::ArrayXd{});
            if (sig.isFinite().all()) {
                signals.row(ii) = sig;
            } else {
                QI::Info("Parameters {} generated non-finite signal", pars.col(ii).transpose());
                signals.row(ii).setZero();
            }
        }
        ref = signals.row(pars.cols() / 2);
    }
    auto const svd = signals.bdcSvd<Eigen::ComputeThinV>();

    QI::Info("Computing projection");
    Eigen::MatrixXd temp  = signals * svd.matrixV().leftCols(B.Get());
    Eigen::MatrixXd proj  = temp * svd.matrixV().leftCols(B.Get()).transpose();
    auto            resid = (signals - proj).stableNorm() / signals.stableNorm();
    QI::Info("Residual {}%", 100 * resid);

    // This will implicitly transpose
    auto bj = json::array();
    for (int ii = 0; ii < B.Get(); ii++) {
        Eigen::VectorXd bc = svd.matrixV().col(ii);
        if (bc.dot(ref) < 0) {
            bc = -bc; // Flip it
        }
        bj.push_back(bc);
    }
    auto sj = json::array();
    for (int ii = 0; ii < signals.rows(); ii++) {
        sj.push_back(signals.row(ii));
    }
    json j = {
        {"basis", bj},
        {"signals", sj},
    };
    QI::WriteJSON(out_path.Get(), j);

    QI::Info("Finished.");
    return EXIT_SUCCESS;
}