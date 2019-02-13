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

    args::HelpFlag       help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag           debug(parser, "DEBUG", "Output debugging messages", {'d', "debug"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "PREFIX", "Add a prefix to output filenames",
                                        {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask",
                                      {'m', "mask"});
    args::ValueFlag<std::string> subregion(
        parser, "REGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK",
        {'s', "subregion"});
    args::Flag resids(parser, "RESIDS", "Write out residuals for each data-point", {'r', "resids"});
    args::ValueFlag<char> algorithm(parser, "ALGO", "Choose algorithm (h)yper/(d)irect, default d",
                                    {'a', "algo"}, 'd');
    args::ValueFlag<std::string> seq_arg(parser, "FILE",
                                         "Read JSON input from file instead of stdin", {"file"});
    args::ValueFlag<float>       simulate(
        parser, "SIMULATE", "Simulate sequence instead of fitting model (argument is noise level)",
        {"simulate"}, 0.0);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::Log(verbose, "Reading sequence information");
    rapidjson::Document json_input = seq_arg ? QI::ReadJSON(seq_arg.Get()) : QI::ReadJSON(std::cin);
    QI::SSFPSequence    ssfp(QI::GetMember(json_input, "SSFP"));
    QI::EllipseModel    model{ssfp};
    if (simulate) {
        QI::SimulateModel<QI::EllipseModel, false>(json_input, model, {}, {QI::CheckPos(ssfp_path)},
                                                   verbose, simulate.Get());
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
        fit_filter->ReadInputs({QI::CheckPos(ssfp_path)}, {}, mask.Get());
        fit_filter->SetBlocks(fit_filter->GetInput(0)->GetNumberOfComponentsPerPixel() /
                              ssfp.size());
        fit_filter->Update();
        fit_filter->WriteOutputs(outarg.Get() + "ES_");
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
