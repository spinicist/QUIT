/*
 *  qidiff.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Args.h"
#include "ImageIO.h"
#include "Util.h"
#include "itkImageRegionConstIterator.h"
#include "itkMultiThreaderBase.h"

//******************************************************************************
// Main
//******************************************************************************
int diff_main(int argc, char **argv) {
    args::ArgumentParser parser("Checks if images are different within a tolerance.\n"
                                "Intended for use with library tests.\n"
                                "http://github.com/spinicist/QUIT");
    args::HelpFlag       help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<std::string> input_path(
        parser, "INPUT", "Input file for difference", {"input"});
    args::ValueFlag<std::string> baseline_path(
        parser, "BASELINE", "Baseline file for difference", {"baseline"});
    args::ValueFlag<double> noise(parser,
                                  "NOISE",
                                  "Added noise level (divide diff by this to get noise factor)",
                                  {"noise"},
                                  0);
    args::Flag              absolute(parser,
                        "ABSOLUTE",
                        "Use absolute difference, not relative (avoids 0/0 problems)",
                        {'a', "abs"});
    QI::ParseArgs(parser, argc, argv, verbose);
    auto input    = QI::ReadImage(QI::CheckValue(input_path), verbose);
    auto baseline = QI::ReadImage(QI::CheckValue(baseline_path), verbose);

    double accum_bias = 0.;
    double accum_var  = 0.;
    auto   mt         = itk::MultiThreaderBase::New();
    mt->SetNumberOfWorkUnits(1);
    mt->ParallelizeImageRegion<3>(
        input->GetBufferedRegion(),
        [&](const QI::VolumeF::RegionType &region) {
            auto input_it = itk::ImageRegionConstIterator<QI::VolumeF>(input, region);
            auto base_it  = itk::ImageRegionConstIterator<QI::VolumeF>(baseline, region);

            for (input_it.GoToBegin(); !input_it.IsAtEnd(); ++input_it, ++base_it) {
                auto const i    = input_it.Get();
                auto const b    = base_it.Get();
                auto const bias = (i - b) / (absolute ? 1. : b);
                accum_bias += bias;
                accum_var += bias * bias;
            }
        },
        nullptr);
    auto const n         = input->GetLargestPossibleRegion().GetNumberOfPixels();
    auto const mean_bias = accum_bias / n;
    auto const var       = accum_var / (n - 1);
    auto const std       = sqrt(var);
    auto const mse       = var + mean_bias * mean_bias;
    auto const me        = sqrt(mse);
    auto const NF        = me / noise.Get();

    if (absolute) {
        QI::Log(verbose,
                "Mean Bias : {}\nVariance  : {}\nSD        : {}\nMSE       : {}\nME        : "
                "{}\nNF           : {}",
                mean_bias,
                var,
                std,
                mse,
                me,
                NF);
    } else {
        QI::Log(verbose,
                "Mean Bias (%): {}\nVariance  (%): {}\nSD        (%): {}\nMSE       (%): {}\nME    "
                "    (%): {}\nNF    "
                "       : {}",
                mean_bias * 100,
                var * 100,
                std * 100,
                mse * 100,
                me * 100,
                NF);
    }
    fmt::print("{}\n", (noise.Get() > 0.) ? NF : me);
    return EXIT_SUCCESS;
}
