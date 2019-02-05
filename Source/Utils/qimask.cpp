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

#include <complex>
#include <string>

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryCrossStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

#include "Args.h"
#include "ImageIO.h"
#include "ImageTypes.h"
#include "Masking.h"
#include "Util.h"

int main(int argc, char **argv) {
    args::ArgumentParser parser("Generates masks in stages.\n"
                                "Stage 1 - Otsu thresholding to generate binary mask\n"
                                "Stage 2 - RATs (optional)\n"
                                "Stage 3 - Hole filling (optional)\n"
                                "http://github.com/spinicist/QUIT");

    args::Positional<std::string> input_path(parser, "INPUT_FILE", "Input file");

    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<std::string> outarg(
        parser, "OUTPUT FILE", "Set output filename, default is input + _mask", {'o', "out"});
    args::ValueFlag<int> volume(
        parser, "VOLUME",
        "Choose volume to mask in multi-volume file. Default 1, -1 selects last volume",
        {'v', "volume"}, 0);
    args::Flag is_complex(parser, "COMPLEX", "Input data is complex, take magnitude first",
                          {'x', "complex"});
    args::ValueFlag<float> lower_threshold(
        parser, "LOWER THRESHOLD",
        "Specify lower intensity threshold for 1st stage, otherwise Otsu's method is used",
        {'l', "lower"}, 0.);
    args::ValueFlag<float> upper_threshold(
        parser, "UPPER THRESHOLD",
        "Specify upper intensity threshold for 1st stage, otherwise Otsu's method is used",
        {'u', "upper"}, std::numeric_limits<float>::infinity());
    args::ValueFlag<float> rats(
        parser, "RATS", "Perform the RATS step, argument is size threshold for connected component",
        {'r', "rats"}, 0.);
    args::ValueFlag<int> fillh_radius(
        parser, "FILL HOLES", "Fill holes in thresholded mask with radius N", {'F', "fillh"}, 0);

    QI::ParseArgs(parser, argc, argv, verbose);

    QI::Log(verbose, "Reading input image: {}", QI::CheckPos(input_path));
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
        QI::Log(verbose, "Masking volume {}...", volume_to_get);
        region.GetModifiableIndex()[3] = volume_to_get;
    } else {
        QI::Fail("Specified mask volume was invalid {} ", volume.Get());
    }
    region.GetModifiableSize()[3] = 0;
    auto vol                      = itk::ExtractImageFilter<QI::SeriesF, QI::VolumeF>::New();
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
        QI::Log(verbose, "Thresholding range: {}-{}", lower_threshold.Get(), upper_threshold.Get());
        mask_image =
            QI::ThresholdMask(intensity_image, lower_threshold.Get(), upper_threshold.Get());
    } else {
        QI::Log(verbose, "Generating Otsu mask");
        mask_image = QI::OtsuMask(intensity_image);
    }

    /*
     *  Stage 2 - RATS
     */
    if (rats) {
        typedef itk::BinaryBallStructuringElement<int, 3>                     TBall;
        typedef itk::BinaryErodeImageFilter<QI::VolumeI, QI::VolumeI, TBall>  TErode;
        typedef itk::BinaryDilateImageFilter<QI::VolumeI, QI::VolumeI, TBall> TDilate;
        float mask_volume  = std::numeric_limits<float>::infinity();
        float voxel_volume = QI::VoxelVolume(mask_image);
        QI::Log(verbose, "Voxel volume: {}", voxel_volume);
        int                  radius = 0;
        QI::VolumeI::Pointer mask_rats;
        while (mask_volume > rats.Get()) {
            radius++;
            TBall           ball;
            TBall::SizeType radii;
            radii.Fill(radius);
            ball.SetRadius(radii);
            ball.CreateStructuringElement();
            TErode::Pointer erode = TErode::New();
            erode->SetInput(mask_image);
            erode->SetForegroundValue(1);
            erode->SetBackgroundValue(0);
            erode->SetKernel(ball);
            erode->Update();
            std::vector<float> kept_sizes = QI::FindLabels(erode->GetOutput(), 0, 1, mask_rats);
            mask_volume                   = kept_sizes[0] * voxel_volume;
            auto dilate                   = TDilate::New();
            dilate->SetKernel(ball);
            dilate->SetInput(mask_rats);
            dilate->SetForegroundValue(1);
            dilate->SetBackgroundValue(0);
            dilate->Update();
            mask_rats = dilate->GetOutput();
            mask_rats->DisconnectPipeline();
            QI::Log(verbose, "Ran RATS iteration, radius = {} volume = {}", radius, mask_volume);
        }
        mask_image = mask_rats;
    }

    /*
     *  Stage 3 - Hole Filling
     */
    QI::VolumeI::Pointer finalMask = ITK_NULLPTR;
    if (fillh_radius) {
        QI::Log(verbose, "Filling holes");
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

    QI::Log(verbose, "Saving mask to: {}", out_path);
    QI::WriteImage(mask_image, out_path);
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
