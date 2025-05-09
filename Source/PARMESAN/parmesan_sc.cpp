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

#define QI_DEBUG_BUILD 1

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
#include "prep_sc.h"

#include <Eigen/Core>
#include <Eigen/SVD>
#include <functional>

int parmesan_fit(args::Subparser &parser) {
    args::Positional<std::string> json_path(parser, "JSON", "Parameter file");
    args::Positional<std::string> input_path(parser, "INPUT", "Input image file");
    args::ValueFlag<std::string>  bpath(parser, "BASIS", "Path to basis JSON file", {'b', "basis"});
    QI_COMMON_ARGS;
    Parse(parser);
    QI::CheckPos(input_path);
    PrepSequence sequence(QI::ReadJSON(json_path.Get())["PrepZTE"]);

    PrepModel model{{}, sequence, ReadBasis(bpath.Get())};
    using FitType = QI::ScaledNumericDiffFit<PrepModel>;
    FitType fit{model};
    auto    fit_filter = QI::ModelFitFilter<FitType>::New(&fit, covar, resids, subregion.Get());
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
    args::ValueFlag<float> noise(parser, "NOISE", "Noise standard deviation", {'n', "noise"}, 0.f);
    args::ValueFlag<std::string> bpath(parser, "BASIS", "Path to basis JSON file", {'b', "basis"});
    QI_CORE_ARGS;
    Parse(parser);
    QI::CheckPos(out_path);
    QI::Info("Reading sequence parameters");
    PrepSequence sequence(QI::ReadJSON(json_path.Get())["PrepZTE"]);
    PrepModel    model{{}, sequence, ReadBasis(bpath.Get())};
    QI::SimulateModel2<decltype(model), false>(
        model, varying_paths.Get(), {}, {out_path.Get()}, mask.Get(), noise.Get(), subregion.Get());
    QI::Info("Finished.");
    return EXIT_SUCCESS;
}

int parmesan_basis(args::Subparser &parser) {
    args::Positional<std::string> json_path(parser, "JSON", "Parameter JSON file");
    args::Positional<std::string> out_path(parser, "OUTPUT", "Basis JSON file");
    args::ValueFlag<int>          B(parser, "B", "Basis size (8)", {'b', "B"}, 8);
    args::ValueFlag<int>          N(parser, "N", "Use N random samples", {'n', "N"}, 16384);

    Parse(parser);
    QI::CheckPos(json_path);
    QI::CheckPos(out_path);
    QI::Info("Reading sequence parameters");
    PrepSequence sequence(QI::ReadJSON(json_path.Get())["PrepZTE"]);
    PrepModel    model{{}, sequence};

    auto const pars =
        N ? QI::RandomPars(model.lo, model.hi, N.Get()) :
            QI::RegularPars(model.lo, model.hi, Eigen::Array<int, 5, 1>{1, 10, 10, 10, 10});

    Eigen::MatrixXd signals(pars.cols(), sequence.size());
    for (int ii = 0; ii < pars.cols(); ii++) {
        signals.row(ii) = model.signal(pars.col(ii), Eigen::ArrayXd{});
    }
    auto const svd = signals.bdcSvd<Eigen::ComputeThinV>();

    QI::Info("Computing projection");
    Eigen::MatrixXd temp  = signals * svd.matrixV().leftCols(B.Get());
    Eigen::MatrixXd proj  = temp * svd.matrixV().leftCols(B.Get()).transpose();
    auto            resid = (signals - proj).stableNorm() / signals.stableNorm();
    QI::Info("Residual {}%", 100 * resid);

    auto bj = json::array();
    for (int ii = 0; ii < B.Get(); ii++) {
        Eigen::VectorXd bc = svd.matrixV().col(ii);
        if (bc.dot(signals.row(pars.cols() / 2)) < 0) {
            bc = -bc; // Flip it
        }
        bj.push_back(bc);
    }
    auto j = json{{"basis", bj}};
    QI::WriteJSON(out_path.Get(), j);

    QI::Info("Finished.");
    return EXIT_SUCCESS;
}

int parmesan_opt(args::Subparser &parser) {
    args::Positional<std::string> json_path(parser, "JSON", "Parameter JSON file");
    // args::Positional<std::string> out_path(parser, "OUTPUT", "Basis JSON file");

    args::ValueFlag<float> tr(parser, "TR", "TR", {"tr"}, 2.e-3);
    args::ValueFlag<float> tramp(parser, "Tramp", "Tramp", {"tramp"}, 10.e-3);
    args::ValueFlag<float> trf(parser, "Trf", "Trf", {"trf"}, 10.e-6);
    args::ValueFlag<int>   sps(parser, "SPS", "SPS", {"sps"}, 64);
    args::ValueFlag<int>   spoils(parser, "SPOILS", "SPOILS", {"spoils"}, 2);
    args::ValueFlag<int>   nP(parser, "NPrep", "Number of preps", {"nprep"}, 8);

    Parse(parser);
    QI::CheckPos(json_path);
    // QI::CheckPos(out_path);
    QI::Info("Reading sequence parameters");

    using ModelType = PrepModel;
    PrepSequence sequence(tr.Get(),
                          tramp.Get(),
                          trf.Get(),
                          0,
                          0,
                          sps.Get(),
                          spoils.Get(),
                          Eigen::ArrayXd::Ones(nP.Get()),
                          Eigen::ArrayXd::Ones(nP.Get()),
                          Eigen::ArrayXd::Ones(nP.Get()) * 1e-3,
                          Eigen::ArrayXd::Zero(nP.Get()));

    ModelType               model{{}, sequence};
    ModelType::VaryingArray v     = (model.lo + model.hi) / 2.;
    v(0)                          = 1; // Set M0 to 1
    ModelType::FixedArray const f = model.fixed_defaults;

    Eigen::MatrixXd FIM(model.NV, model.NV);
    for (int ij = 0; ij < model.NV; ij++) {
        auto const dj = model.dsdθ(v, f, ij);
        for (int ii = 0; ii < model.NV; ii++) {
            auto const di = model.dsdθ(v, f, ii);
            FIM(ii, ij)   = di.matrix().dot(dj.matrix());
        }
    }

    QI_DBMAT(FIM);
    Eigen::MatrixXd const iFIM = FIM.inverse();
    Eigen::VectorXd const CRB  = iFIM.diagonal();
    Eigen::VectorXd const nCRB = CRB.array() * sequence.TR * sequence.SPS / v.square();
    QI_DBMAT(iFIM);
    QI_DBVEC(CRB);
    QI_DBVEC(nCRB);
    // QI_DBVEC(cov);

    QI::Info("Finished.");
    return EXIT_SUCCESS;
}
