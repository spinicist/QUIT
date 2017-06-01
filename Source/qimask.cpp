/*
 *  qimask.cpp
 *
 *  Created by Tobias Wood on 2015/09/22.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <string>
#include <complex>

#include "itkComplexToModulusImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

#include "QI/Args.h"
#include "QI/Types.h"
#include "QI/Util.h"
#include "QI/IO.h"

int main(int argc, char **argv) {
    args::ArgumentParser parser(
    "Generates masks in stages.\n"
    "Stage 1 - Otsu thresholding to generate binary mask\n"
    "Stage 2 - Fill holes (optional)\n"
    "http://github.com/spinicist/QUIT"
    );

    args::Positional<std::string> input_path(parser, "INPUT_FILE", "Input file");

    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<std::string> outarg(parser, "OUTPUT FILE", "Set output filename, default is input + _mask", {'o', "out"});
    args::ValueFlag<int> volume(parser, "VOLUME", "Choose volume to mask in multi-volume file. Default 1, -1 selects last volume", {'v', "volume"}, 0);
    args::Flag     is_complex(parser, "COMPLEX", "Input data is complex, take magnitude first", {'x', "complex"});
    args::ValueFlag<int> fillh_radius(parser, "FILL HOLES", "Fill holes in thresholded mask with radius N", {'F', "fillh"}, 0);
    args::Flag     run_bet(parser, "RUN BET", "Run the Brain Extraction Tool stage", {'B', "bet"});

    QI::ParseArgs(parser, argc, argv);

    if (verbose) std::cout << "Reading input image: " << QI::CheckPos(input_path) << std::endl;
    QI::SeriesF::Pointer vols = ITK_NULLPTR;
    if (is_complex) {
        vols = QI::ReadMagnitudeImage<QI::SeriesF>(QI::CheckPos(input_path));
    } else {
        vols = QI::ReadImage<QI::SeriesF>(QI::CheckPos(input_path));
    }

    std::string out_path;
    std::cout << "outarg is: " << outarg.Get() << std::endl;
    if (outarg.Get() == "") {
        out_path = QI::StripExt(input_path.Get()) + "_mask" + QI::OutExt();
    } else {
        out_path = outarg.Get();
    }

    // Extract one volume to process
    auto region = vols->GetLargestPossibleRegion();
    if (std::abs(volume.Get()) < region.GetSize()[3]) {
        region.GetModifiableIndex()[3] = (region.GetSize()[3] + volume.Get()) % region.GetSize()[3];
        if (verbose) std::cout << "Using volume " << region.GetIndex()[3] << " out of " << region.GetSize()[3] << std::endl;
    } else {
        std::cerr << "Specified volume was invalid: " << volume.Get() << std::endl;
        return EXIT_FAILURE;
    }
    region.GetModifiableSize()[3] = 0;
    auto vol = itk::ExtractImageFilter<QI::SeriesF, QI::VolumeF>::New();
    vol->SetExtractionRegion(region);
    vol->SetInput(vols);
    vol->SetDirectionCollapseToSubmatrix();

    /*
     *  Stage 1 - Otsu or Threshold
     */
    // Use an Otsu Threshold filter to generate the mask
    if (verbose) std::cout << "Generating Otsu mask" << std::endl;
    auto otsuFilter = itk::OtsuThresholdImageFilter<QI::VolumeF, QI::VolumeF>::New();
    otsuFilter->SetInput(vol->GetOutput());
    otsuFilter->SetOutsideValue(1);
    otsuFilter->SetInsideValue(0);
    otsuFilter->Update();
    
    /*
     *  Stage 2 - Hole Filling
     */
    QI::VolumeF::Pointer finalMask = ITK_NULLPTR;
    if (fillh_radius.Get() > 0) {
        if (verbose) std::cout << "Filling holes" << std::endl;
        auto fillHoles = itk::VotingBinaryIterativeHoleFillingImageFilter<QI::VolumeF>::New();
        itk::VotingBinaryIterativeHoleFillingImageFilter<QI::VolumeF>::InputSizeType radius;
        radius.Fill(fillh_radius.Get());
        fillHoles->SetInput(otsuFilter->GetOutput());
        fillHoles->SetRadius(radius);
        fillHoles->SetMajorityThreshold(2); // Corresponds to (rad^3-1)/2 + 2 threshold
        fillHoles->SetBackgroundValue(0.0);
        fillHoles->SetForegroundValue(1.0);
        fillHoles->SetMaximumNumberOfIterations(3);
        fillHoles->Update();
        finalMask = fillHoles->GetOutput();
        finalMask->DisconnectPipeline();
    } else {
        finalMask = otsuFilter->GetOutput();
    }
    if (verbose) std::cout << "Saving mask to: " << out_path << std::endl;
    QI::WriteImage(finalMask, out_path);
    if (verbose) std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}

