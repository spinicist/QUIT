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
#include "itkBinaryCrossStructuringElement.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkExtractImageFilter.h"

#include "Args.h"
#include "ImageTypes.h"
#include "Util.h"
#include "ImageIO.h"
#include "Masking.h"

int main(int argc, char **argv) {
    args::ArgumentParser parser(
    "Generates masks in stages.\n"
    "Stage 1 - Otsu thresholding to generate binary mask\n"
    "Stage 2 - RATs (optional)\n"
    "Stage 3 - Hole filling (optional)\n"
    "http://github.com/spinicist/QUIT"
    );

    args::Positional<std::string> input_path(parser, "INPUT_FILE", "Input file");

    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<std::string> outarg(parser, "OUTPUT FILE", "Set output filename, default is input + _mask", {'o', "out"});
    args::ValueFlag<int> volume(parser, "VOLUME", "Choose volume to mask in multi-volume file. Default 1, -1 selects last volume", {'v', "volume"}, 0);
    args::Flag     is_complex(parser, "COMPLEX", "Input data is complex, take magnitude first", {'x', "complex"});
    args::ValueFlag<float> lower_threshold(parser, "LOWER THRESHOLD", "Specify lower intensity threshold for 1st stage, otherwise Otsu's method is used", {'l', "lower"}, 0.);
    args::ValueFlag<float> upper_threshold(parser, "UPPER THRESHOLD", "Specify upper intensity threshold for 1st stage, otherwise Otsu's method is used", {'u', "upper"}, std::numeric_limits<float>::infinity());
    args::ValueFlag<float> rats(parser, "RATS", "Perform the RATS step, argument is size threshold for connected component", {'r', "rats"}, 0.);
    args::ValueFlag<int> fillh_radius(parser, "FILL HOLES", "Fill holes in thresholded mask with radius N", {'F', "fillh"}, 0);

    QI::ParseArgs(parser, argc, argv, verbose);

    if (verbose) std::cout << "Reading input image: " << QI::CheckPos(input_path) << std::endl;
    QI::SeriesF::Pointer vols = ITK_NULLPTR;
    if (is_complex) {
        vols = QI::ReadMagnitudeImage<QI::SeriesF>(QI::CheckPos(input_path));
    } else {
        vols = QI::ReadImage<QI::SeriesF>(QI::CheckPos(input_path));
    }

    std::string out_path;
    if (outarg.Get() == "") {
        out_path = QI::StripExt(input_path.Get()) + "_mask" + QI::OutExt();
    } else {
        out_path = outarg.Get();
    }

    // Extract one volume to process
    auto region = vols->GetLargestPossibleRegion();
    if (static_cast<size_t>(std::abs(volume.Get())) < region.GetSize()[3]) {
        int volume_to_get = (region.GetSize()[3] + volume.Get()) % region.GetSize()[3];
        if (verbose) std::cout << "Using volume " << volume_to_get << std::endl;
        region.GetModifiableIndex()[3] = volume_to_get;
    } else {
        std::cerr << "Specified volume was invalid: " << volume.Get() << std::endl;
        return EXIT_FAILURE;
    }
    region.GetModifiableSize()[3] = 0;
    auto vol = itk::ExtractImageFilter<QI::SeriesF, QI::VolumeF>::New();
    vol->SetExtractionRegion(region);
    vol->SetInput(vols);
    vol->SetDirectionCollapseToSubmatrix();
    vol->Update();
    QI::VolumeF::Pointer intensity_image = vol->GetOutput();
    intensity_image->DisconnectPipeline();

    /*
     *  Stage 1 - Otsu or Threshold
     */
    QI::VolumeI::Pointer mask_image = ITK_NULLPTR;
    if (lower_threshold || upper_threshold) {
        if (verbose) std::cout << "Thresholding range: " << lower_threshold.Get() << "-" << upper_threshold.Get() << std::endl;
        mask_image = QI::ThresholdMask(intensity_image, lower_threshold.Get(), upper_threshold.Get());
    } else {
        if (verbose) std::cout << "Generating Otsu mask" << std::endl;
        mask_image = QI::OtsuMask(intensity_image);
    }

    /*
     *  Stage 2 - RATS
     */
     if (rats) {
        typedef itk::BinaryBallStructuringElement<int, 3> TBall;
        typedef itk::BinaryErodeImageFilter<QI::VolumeI, QI::VolumeI, TBall> TErode;
        typedef itk::BinaryDilateImageFilter<QI::VolumeI, QI::VolumeI, TBall> TDilate;
        float mask_volume = std::numeric_limits<float>::infinity();
        float voxel_volume = QI::VoxelVolume(mask_image);
        if (verbose) std::cout << "Voxel volume: " << voxel_volume << std::endl;
        int radius = 0;
        QI::VolumeI::Pointer mask_rats;
        while (mask_volume > rats.Get()) {
            radius++;
            TBall ball;
            TBall::SizeType radii; radii.Fill(radius);
            ball.SetRadius(radii);
            ball.CreateStructuringElement();
            TErode::Pointer erode = TErode::New();
            erode->SetInput(mask_image);
            erode->SetForegroundValue(1);
            erode->SetBackgroundValue(0);
            erode->SetKernel(ball);
            erode->Update();
            std::vector<float> kept_sizes = QI::FindLabels(erode->GetOutput(), 0, 1, mask_rats);
            mask_volume = kept_sizes[0] * voxel_volume;
            auto dilate = TDilate::New();
            dilate->SetKernel(ball);
            dilate->SetInput(mask_rats);
            dilate->SetForegroundValue(1);
            dilate->SetBackgroundValue(0);
            dilate->Update();
            mask_rats = dilate->GetOutput();
            mask_rats->DisconnectPipeline();
            if (verbose) std::cout << "Ran RATS iteration, current radius: " << radius << " volume is: " << mask_volume << std::endl;
         }
         mask_image = mask_rats;
     }

    /*
     *  Stage 3 - Hole Filling
     */
    QI::VolumeI::Pointer finalMask = ITK_NULLPTR;
    if (fillh_radius) {
        if (verbose) std::cout << "Filling holes" << std::endl;
        auto fillHoles = itk::VotingBinaryIterativeHoleFillingImageFilter<QI::VolumeI>::New();
        itk::VotingBinaryIterativeHoleFillingImageFilter<QI::VolumeI>::InputSizeType radius;
        radius.Fill(fillh_radius.Get());
        fillHoles->SetInput(mask_image);
        fillHoles->SetRadius(radius);
        fillHoles->SetMajorityThreshold(2); // Corresponds to (rad^3-1)/2 + 2 threshold
        fillHoles->SetBackgroundValue(0);
        fillHoles->SetForegroundValue(1);
        fillHoles->SetMaximumNumberOfIterations(3);
        fillHoles->Update();
        mask_image = fillHoles->GetOutput();
        mask_image->DisconnectPipeline();
    }

    if (verbose) std::cout << "Saving mask to: " << out_path << std::endl;
    QI::WriteImage(mask_image, out_path);
    if (verbose) std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}

