/*
 *  qi_ellipse.cpp
 *
 *  Copyright (c) 2016, 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>

#include "Args.h"
#include "DirectFit.h"
#include "HyperFit.h"
#include "ImageIO.h"
#include "ModelFitFilter.h"
#include "SimulateModel.h"
#include "Util.h"

int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser(
        "Calculates the ellipse parameters G,a,b,f0 & psi0 from SSFP data.\nInput must be a single "
        "complex image with at least 6 phase increments.\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> ssfp_path(parser, "SSFP_FILE", "Input SSFP file");

    QI_COMMON_ARGS;
    args::Flag            debug(parser, "DEBUG", "Output debugging messages", {'d', "debug"});
    args::ValueFlag<char> algorithm(
        parser, "ALGO", "Choose algorithm (h)yper/(d)irect, default d", {'a', "algo"}, 'd');
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::CheckPos(ssfp_path);
    QI::Log(verbose, "Reading sequence information");
    json input = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    auto ssfp  = input.at("SSFP").get<QI::SSFPSequence>();
    QI::Log(verbose, "{}", ssfp);
    QI::EllipseModel model{ssfp};
    if (simulate) {
        QI::SimulateModel<QI::EllipseModel, false>(
            input, model, {}, {ssfp_path.Get()}, verbose, simulate.Get());
    } else {
        QI::EllipseFit *fit = nullptr;
        switch (algorithm.Get()) {
        case 'h':
            fit = new QI::HyperFit(model);
            break;
        case 'd':
            fit = new QI::DirectFit(model);
            break;
        }
        auto fit_filter =
            QI::ModelFitFilter<QI::EllipseFit>::New(fit, verbose, resids, subregion.Get());
        fit_filter->ReadInputs({ssfp_path.Get()}, {}, mask.Get());
        fit_filter->SetBlocks(fit_filter->GetInput(0)->GetNumberOfComponentsPerPixel() /
                              ssfp.size());
        fit_filter->Update();
        fit_filter->WriteOutputs(prefix.Get() + "ES_");
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
