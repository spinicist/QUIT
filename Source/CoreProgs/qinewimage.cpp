/*
 *  qinewimage.cpp
 *
 *  Created by Tobias Wood on 02/06/2015.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <random>

#include "itkImage.h"
#include "itkImageLinearIteratorWithIndex.h"

#include "Args.h"
#include "ImageIO.h"
#include "Util.h"

template <int dim>
void make_image(args::Positional<std::string> &fName,
                args::ValueFlag<std::string> & size_arg,
                args::ValueFlag<std::string> & spacing_arg,
                args::ValueFlag<std::string> & origin_arg,
                args::ValueFlag<float> &       fill_arg,
                args::ValueFlag<int> &         grad_dim,
                args::ValueFlag<std::string> & grad_arg,
                args::ValueFlag<int> &         steps_arg,
                args::ValueFlag<float> &       wrap_arg) {
    typedef itk::Image<float, dim> ImageType;
    using RegionType  = typename ImageType::RegionType;
    using IndexType   = typename ImageType::IndexType;
    using SizeType    = typename ImageType::SizeType;
    using SpacingType = typename ImageType::SpacingType;
    using PointType   = typename ImageType::PointType;
    auto       newimg = ImageType::New();
    RegionType imgRegion;
    IndexType  imgIndex;
    imgIndex.Fill(0);
    if (grad_dim.Get() >= dim)
        QI::Fail("Fill dimension is larger than image dimension");
    SizeType imgSize;
    if (size_arg) {
        QI::ArrayArg<SizeType, dim>(size_arg.Get(), imgSize);
    } else {
        imgSize.Fill(1);
    }
    SpacingType imgSpacing;
    if (spacing_arg) {
        QI::ArrayArg<SpacingType, dim>(spacing_arg.Get(), imgSpacing);
    } else {
        imgSpacing.Fill(1.0);
    }
    PointType imgOrigin;
    if (origin_arg) {
        QI::ArrayArg<PointType, dim>(origin_arg.Get(), imgOrigin);
    } else {
        for (int i = 0; i < 3; i++) {
            imgOrigin[i] = -imgSpacing[i] * (imgSize[i] - 1) / 2.0;
        }
    }
    QI::Log(verbose,
            "Dimensions: {} Size: {} Spacing: {} Origin: {}",
            dim,
            imgSize,
            imgSpacing,
            imgOrigin);
    Eigen::Array2f end_vals{0, 0};

    float deltaVal    = 0;
    int   step_length = 1;
    if (grad_arg) {
        QI::ArrayArgF<Eigen::Array2f, 2>(grad_arg.Get(), end_vals);
        QI::Log(verbose, "Fill starts at {}, ends at {}", end_vals[0], end_vals[1]);
        if (steps_arg) {
            deltaVal    = (end_vals[1] - end_vals[0]) / steps_arg.Get();
            step_length = (imgSize[grad_dim.Get()] - 1) / steps_arg.Get();
            QI::Log(verbose,
                    "Number of steps is {}, step length is {}, step value is {}",
                    steps_arg.Get(),
                    step_length,
                    deltaVal);
        } else {
            deltaVal = (end_vals[1] - end_vals[0]) / (imgSize[grad_dim.Get()] - 1);
            QI::Log(verbose, "Smooth gradient, delta value {}", deltaVal);
        }
        QI::Log(verbose, "Along dimension {}", grad_dim.Get());

    } else if (fill_arg) {
        end_vals = fill_arg.Get();
        deltaVal = 0;
        QI::Log(verbose, "Fill with constant value: {}", fill_arg.Get());
    }

    imgRegion.SetIndex(imgIndex);
    imgRegion.SetSize(imgSize);
    newimg->SetRegions(imgRegion);
    newimg->SetSpacing(imgSpacing);
    newimg->SetOrigin(imgOrigin);
    newimg->Allocate();

    itk::ImageLinearIteratorWithIndex<ImageType> it(newimg, imgRegion);
    it.SetDirection(grad_dim.Get());
    QI::Log(verbose, "Filling...");
    it.GoToBegin();
    while (!it.IsAtEnd()) {
        float val = end_vals[0];
        while (!it.IsAtEndOfLine()) {
            it.Set(val);
            ++it;
            if ((it.GetIndex()[grad_dim.Get()] % step_length) == 0) {
                val += deltaVal;
            }
            if (wrap_arg) {
                val = fmod(val, wrap_arg.Get());
            }
        }
        it.NextLine();
    }
    QI::Log(verbose, "Writing file to: {}", QI::CheckPos(fName));
    QI::WriteImage(newimg, QI::CheckPos(fName), verbose);
}

int newimage_main(args::Subparser &parser) {
    args::Positional<std::string> fName(parser, "OUTPUT", "Output filename");

    args::ValueFlag<int> dims_arg(parser, "DIMS", "Image dimension, default 3", {'d', "dims"}, 3);
    args::ValueFlag<std::string> size_arg(parser, "SIZE", "Image size", {'s', "size"});
    args::ValueFlag<std::string> spacing_arg(parser, "SPACING", "Voxel spacing", {'p', "spacing"});
    args::ValueFlag<std::string> origin_arg(parser, "ORIGIN", "Image origin", {'o', "origin"});
    args::ValueFlag<float>       fill_arg(parser, "FILL", "Fill with value", {'f', "fill"});
    args::ValueFlag<int>         grad_dim(
        parser, "GRAD DIM", "Fill with a gradient along this axis", {'g', "grad_dim"}, 0);
    args::ValueFlag<std::string> grad_arg(
        parser, "GRAD END", "Gradient low/high values (low, high)", {'v', "grad_vals"}, "0,0");
    args::ValueFlag<int> steps_arg(
        parser, "STEPS", "Number of discrete steps (steps)", {'t', "steps"}, 1);
    args::ValueFlag<float> wrap_arg(parser, "WRAP", "Wrap image values", {'w', "wrap"});
    parser.Parse();

    if (dims_arg.Get() == 3) {
        make_image<3>(fName,
                      size_arg,
                      spacing_arg,
                      origin_arg,
                      fill_arg,
                      grad_dim,
                      grad_arg,
                      steps_arg,
                      wrap_arg);
    } else if (dims_arg.Get() == 4) {
        make_image<4>(fName,
                      size_arg,
                      spacing_arg,
                      origin_arg,
                      fill_arg,
                      grad_dim,
                      grad_arg,
                      steps_arg,
                      wrap_arg);
    } else {
        std::cerr << "Unsupported dimension: " << dims_arg.Get() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
