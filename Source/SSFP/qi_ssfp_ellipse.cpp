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
    args::ValueFlag<char> algorithm(parser, "ALGO", "Choose algorithm (h)yper/(d)irect, default d",
                                    {'a', "algo"}, 'd');
    args::ValueFlag<std::string> seq_arg(parser, "FILE",
                                         "Read JSON input from file instead of stdin", {"file"});
    args::ValueFlag<float>       simulate(
        parser, "SIMULATE", "Simulate sequence instead of fitting model (argument is noise level)",
        {"simulate"}, 0.0);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::CheckPos(ssfp_path);
    QI::Log(verbose, "Reading sequence information");
    rapidjson::Document json_input = seq_arg ? QI::ReadJSON(seq_arg.Get()) : QI::ReadJSON(std::cin);
    QI::SSFPSequence    ssfp(QI::GetMember(json_input, "SSFP"));
    QI::EllipseModel    model{ssfp};
    if (simulate) {
        QI::SimulateModel<QI::EllipseModel, false>(json_input, model, {}, {ssfp_path.Get()},
                                                   verbose, simulate.Get());
    } else {
        auto input = QI::ReadVectorImage<std::complex<float>>(QI::CheckPos(ssfp_path), verbose);
        QI::EllipseFit *fit = nullptr;
        switch (algorithm.Get()) {
        case 'h':
            fit = new QI::HyperFit(model);
            break;
        case 'd':
            fit = new QI::DirectFit(model);
            break;
        }
        auto fit_filter = QI::ModelFitFilter<QI::EllipseFit>::New(fit, verbose, false);
        fit_filter->SetInput(0, input);
        fit_filter->SetBlocks(input->GetNumberOfComponentsPerPixel() / ssfp.size());
        if (mask)
            fit_filter->SetMask(QI::ReadImage(mask.Get(), verbose));
        if (subregion)
            fit_filter->SetSubregion(QI::RegionArg(args::get(subregion)));
        fit_filter->Update();
        std::string outPrefix = outarg.Get() + "ES_";
        for (int i = 0; i < QI::EllipseModel::NV; i++) {
            QI::WriteVectorImage(fit_filter->GetOutput(i),
                                 outPrefix + QI::EllipseModel::varying_names.at(i) + QI::OutExt());
        }
        // QI::WriteImage(fit_filter->GetResidualOutput(), outPrefix + "residual" + QI::OutExt());
        // if (resids) {
        //     QI::WriteVectorImage(fit_filter->GetResidualsOutput(0), outPrefix + "all_residuals" +
        //     QI::OutExt());
        // }
        // if (its) {
        //     QI::WriteImage(fit_filter->GetFlagOutput(), outPrefix + "iterations" + QI::OutExt());
        // }
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
