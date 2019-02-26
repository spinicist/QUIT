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
#include "itkDivideImageFilter.h"
#include "itkSquareImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkSubtractImageFilter.h"

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
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

    auto diffFilter = itk::SubtractImageFilter<QI::VolumeF>::New();
    diffFilter->SetInput1(input);
    diffFilter->SetInput2(baseline);

    auto sqr_norm = itk::SquareImageFilter<QI::VolumeF, QI::VolumeF>::New();
    if (absolute) {
        sqr_norm->SetInput(diffFilter->GetOutput());
    } else {
        auto diff_norm = itk::DivideImageFilter<QI::VolumeF, QI::VolumeF, QI::VolumeF>::New();
        diff_norm->SetInput1(diffFilter->GetOutput());
        diff_norm->SetInput2(baseline);
        diff_norm->Update();
        sqr_norm->SetInput(diff_norm->GetOutput());
    }
    auto stats = itk::StatisticsImageFilter<QI::VolumeF>::New();
    stats->SetInput(sqr_norm->GetOutput());
    stats->Update();

    const double mean_sqr_diff      = stats->GetMean();
    const double root_mean_sqr_diff = sqrt(mean_sqr_diff);
    const double diff = (noise.Get() > 0) ? root_mean_sqr_diff / noise.Get() : root_mean_sqr_diff;
    QI::Log(verbose,
            "Mean Square Diff: {}\nSquare-root mean square diff: {}\nRelative noise: {}\nRelative "
            "Diff: {}",
            mean_sqr_diff,
            root_mean_sqr_diff,
            noise.Get(),
            diff);
    fmt::print("{}", diff);
}
