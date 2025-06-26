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
#include "FitScaledNumericMultiStart.h"
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
    args::Positional<std::string> jpath(parser, "JSON", "Parameter file");
    args::Positional<std::string> input_path(parser, "INPUT", "Input image file");
    args::ValueFlag<std::string>  bpath(parser, "BASIS", "Path to basis JSON file", {'b', "basis"});
    args::Flag                    nodf(parser, "F", "No off-resonance", {'f', "nodf"});
    args::Flag                    gauss(parser, "G", "Gauss Prep Pulses", {'g', "gauss"});
    QI_COMMON_ARGS;
    Parse(parser);
    QI::CheckPos(input_path);
    PrepSequence sequence(QI::ReadJSON(jpath.Get())["PrepZTE"]);
    if (nodf) {
        PrepModel2 model{{}, sequence, ReadBasis(bpath.Get())};
        using FitType = QI::ScaledNumericDiffFit<PrepModel2>;
        FitType fit{model};
        auto    fit_filter = QI::ModelFitFilter<FitType>::New(&fit, covar, resids, subregion.Get());
        fit_filter->ReadInputs({input_path.Get()}, {}, mask.Get());
        fit_filter->SetRequestedRegionFromString(subregion.Get());
        fit_filter->Update();
        fit_filter->WriteOutputs(prefix.Get() + "PARMESAN_");
    } else {
        PrepModel model{{}, sequence, ReadBasis(bpath.Get()), gauss};
        using FitType = QI::ScaledNumericDiffMultiStartFit<PrepModel>;
        Eigen::ArrayXd f0_starts(1);
        f0_starts << 0.;
        FitType fit{model, f0_starts};
        auto    fit_filter = QI::ModelFitFilter<FitType>::New(&fit, covar, resids, subregion.Get());
        fit_filter->ReadInputs({input_path.Get()}, {}, mask.Get());
        fit_filter->SetRequestedRegionFromString(subregion.Get());
        fit_filter->Update();
        fit_filter->WriteOutputs(prefix.Get() + "PARMESAN_");
    }
    QI::Info("Finished.");
    return EXIT_SUCCESS;
}

int parmesan_sim(args::Subparser &parser) {
    args::Positional<std::string>     jpath(parser, "JSON", "Parameter file");
    args::Positional<std::string>     opath(parser, "OUTPUT", "Simulation output file");
    args::PositionalList<std::string> vpaths(parser, "INPUT", "Input parameter maps");
    args::ValueFlag<float> noise(parser, "NOISE", "Noise standard deviation", {'n', "noise"}, 0.f);
    args::ValueFlag<std::string> bpath(parser, "BASIS", "Path to basis JSON file", {'b', "basis"});
    args::Flag                   nodf(parser, "F", "No off-resonance", {'f', "nodf"});
    args::Flag                   gauss(parser, "G", "Gauss Prep Pulses", {'g', "gauss"});
    QI_CORE_ARGS;
    Parse(parser);
    QI::CheckPos(opath);
    QI::Info("Reading sequence parameters");
    PrepSequence sequence(QI::ReadJSON(jpath.Get())["PrepZTE"]);
    if (nodf) {
        PrepModel2 model{{}, sequence, ReadBasis(bpath.Get())};
        QI::SimulateModel2<decltype(model), false>(
            model, vpaths.Get(), {}, {opath.Get()}, mask.Get(), noise.Get(), subregion.Get());
    } else {
        PrepModel model{{}, sequence, ReadBasis(bpath.Get()), gauss};
        QI::SimulateModel2<decltype(model), false>(
            model, vpaths.Get(), {}, {opath.Get()}, mask.Get(), noise.Get(), subregion.Get());
    }
    QI::Info("Finished.");
    return EXIT_SUCCESS;
}

template <typename ModelT>
auto RunBasis(PrepSequence const &sequence, int const N, bool const gauss) -> Eigen::MatrixXd {
    ModelT model{{}, sequence, Eigen::MatrixXd(), gauss};

    auto const pars =
        (N > 0) ? QI::RandomPars(model.lo, model.hi, N) :
                  QI::RegularPars(model.lo, model.hi, Eigen::Array<int, 5, 1>{1, 10, 10, 10, 10});

    Eigen::MatrixXd signals(pars.cols(), sequence.size());
    for (int ii = 0; ii < pars.cols(); ii++) {
        signals.row(ii) = model.signal(pars.col(ii), Eigen::ArrayXd{});
    }
    return signals;
}

int parmesan_basis(args::Subparser &parser) {
    args::Positional<std::string> jpath(parser, "JSON", "Parameter JSON file");
    args::Positional<std::string> opath(parser, "OUTPUT", "Basis JSON file");
    args::ValueFlag<int>          B(parser, "B", "Basis size (8)", {'b', "B"}, 8);
    args::ValueFlag<int>          N(parser, "N", "Use N random samples", {'n', "N"}, 0);
    args::Flag                    nodf(parser, "F", "No off-resonance", {'f', "nodf"});
    args::Flag                    gauss(parser, "G", "Gauss Prep Pulses", {'g', "gauss"});
    Parse(parser);
    QI::CheckPos(jpath);
    QI::CheckPos(opath);

    PrepSequence sequence(QI::ReadJSON(jpath.Get())["PrepZTE"]);
    auto const   signals = nodf ? RunBasis<PrepModel2>(sequence, N.Get(), gauss) :
                                  RunBasis<PrepModel>(sequence, N.Get(), gauss);
    auto const   svd     = signals.bdcSvd<Eigen::ComputeThinV>();

    QI::Info("Computing projection");
    Eigen::MatrixXd temp  = signals * svd.matrixV().leftCols(B.Get());
    Eigen::MatrixXd proj  = temp * svd.matrixV().leftCols(B.Get()).transpose();
    auto            resid = (signals - proj).stableNorm() / signals.stableNorm();
    QI::Info("Residual {}%", 100 * resid);

    auto bj = json::array();
    for (int ii = 0; ii < B.Get(); ii++) {
        Eigen::VectorXd bc = svd.matrixV().col(ii);
        if (bc.dot(signals.row(signals.rows() / 2)) < 0) {
            bc = -bc; // Flip it
        }
        bj.push_back(bc);
    }
    auto j = json{{"basis", bj}};
    QI::WriteJSON(opath.Get(), j);

    QI::Info("Finished.");
    return EXIT_SUCCESS;
}
