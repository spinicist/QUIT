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

#include <iostream>

#include "itkSubtractImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkSquareImageFilter.h"
#include "Util.h"
#include "Args.h"
#include "ImageIO.h"

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    args::ArgumentParser parser("Checks if images are different within a tolerance.\n"
                                "Intended for use with library tests.\n"
                                "http://github.com/spinicist/QUIT");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<std::string> input_path(parser, "INPUT", "Input file for difference", {"input"});
    args::ValueFlag<std::string> baseline_path(parser, "BASELINE", "Baseline file for difference", {"baseline"});
    args::ValueFlag<double>      tolerance(parser, "TOLERANCE", "Tolerance (mean percent difference)", {"tolerance"}, 0);
    args::ValueFlag<double>      noise(parser, "NOISE", "Added noise level, tolerance is relative to this", {"noise"}, 1);
    args::Flag                   absolute(parser, "ABSOLUTE", "Use absolute difference, not relative (avoids 0/0 problems)", {'a', "abs"});
    QI::ParseArgs(parser, argc, argv, verbose);
    if (verbose) std::cout << "Reading input: " << QI::CheckValue(input_path) << std::endl;
    if (verbose) std::cout << "Reading baseline: " << QI::CheckValue(baseline_path) << std::endl;
    auto input = QI::ReadImage(QI::CheckValue(input_path));
    auto baseline = QI::ReadImage(QI::CheckValue(baseline_path));
    
    auto diff = itk::SubtractImageFilter<QI::VolumeF>::New();
    diff->SetInput1(input);
    diff->SetInput2(baseline);

    auto sqr_norm = itk::SquareImageFilter<QI::VolumeF, QI::VolumeF>::New();
    if (absolute) {
        sqr_norm->SetInput(diff->GetOutput());
    } else {
        auto diff_norm = itk::DivideImageFilter<QI::VolumeF, QI::VolumeF, QI::VolumeF>::New();
        diff_norm->SetInput1(diff->GetOutput());
        diff_norm->SetInput2(baseline);
        diff_norm->Update();
        sqr_norm->SetInput(diff_norm->GetOutput());
    }
    auto stats = itk::StatisticsImageFilter<QI::VolumeF>::New();
    stats->SetInput(sqr_norm->GetOutput());
    stats->Update();

    const double mean_sqr_diff = stats->GetMean();
    const double root_mean_sqr_diff = sqrt(mean_sqr_diff);
    const double rel_diff = root_mean_sqr_diff / noise.Get();
    const bool passed = rel_diff <= tolerance.Get();
    if (verbose) {
        std::cout << "Mean Square Diff: " << mean_sqr_diff
                  << "\nRelative noise: " << noise.Get()
                  << "\nSquare-root mean square diff: " << root_mean_sqr_diff
                  << "\nRelative Diff: " << rel_diff
                  << "\nTolerance: " << tolerance.Get()
                  << "\nResult: " << (passed ? "Passed" : "Failed") << std::endl;
    }
    if (passed) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }
}


