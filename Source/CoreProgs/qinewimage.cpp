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

/*
 * Define arguments globally because I am lazy
 */
args::ArgumentParser
    parser("This is a tool to create 3D Nifti files, either blank headers with orientation "
           "information, e.g. for registration, or files filled with simple patterns of "
           "data e.g. solid values, gradients, or blocks. The default is to create a 3D file "
           "filled with zeros. Choose from the options below to create something else."
           "\nhttp://github.com/spinicist/QUIT");

args::Positional<std::string> fName(parser, "OUTPUT", "Output filename");

args::HelpFlag       help(parser, "HELP", "Show this help menu", {'h', "help"});
args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
args::ValueFlag<int> dims_arg(parser, "DIMS", "Image dimension, default 3", {'d', "dims"}, 3);
args::ValueFlag<std::string> size_arg(parser, "SIZE", "Image size", {'s', "size"});
args::ValueFlag<std::string> spacing_arg(parser, "SPACING", "Voxel spacing", {'p', "spacing"});
args::ValueFlag<std::string> origin_arg(parser, "ORIGIN", "Image origin", {'o', "origin"});
args::ValueFlag<std::string> fill_arg(parser, "FILL", "Fill with value", {'f', "fill"});
args::ValueFlag<std::string> grad_arg(parser, "GRADIENT", "Fill with gradient (dim, low, high)",
                                      {'g', "grad"});
args::ValueFlag<std::string>
                        steps_arg(parser, "STEPS", "Fill with discrete steps (dim, low, high, steps)", {'t', "step"});
args::ValueFlag<double> wrap_arg(parser, "WRAP", "Wrap image values", {'w', "wrap"});

template <int dim> void make_image() {
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
    SizeType imgSize;
    imgSize.Fill(1);
    SpacingType imgSpacing;
    imgSpacing.Fill(1.0);
    PointType imgOrigin;
    imgOrigin.Fill(0.);

    if (size_arg)
        QI::ArrayArg<SizeType, dim>(size_arg.Get(), imgSize);
    if (spacing_arg)
        QI::ArrayArg<SpacingType, dim>(spacing_arg.Get(), imgSpacing);
    if (origin_arg)
        QI::ArrayArg<PointType, dim>(origin_arg.Get(), imgOrigin);
    QI::Log(verbose, "Size: {} Spacing: {} Origin: {}", imgSize, imgSpacing, imgOrigin);
    enum class FillTypes { Fill, Gradient, Steps };
    FillTypes fillType = FillTypes::Fill;
    float     startVal = 0, deltaVal = 0, stopVal = 0;
    int       stepLength = 1, fillDim = 0;
    if (fill_arg) {
        fillType = FillTypes::Fill;
        startVal = std::stod(fill_arg.Get());
        QI::Log(verbose, "Fill with constant value: {}", startVal);
    } else if (grad_arg) {
        fillType = FillTypes::Gradient;
        std::istringstream stream(grad_arg.Get());
        stream >> fillDim;
        stream >> startVal;
        stream >> stopVal;
        deltaVal   = (stopVal - startVal) / (imgSize[fillDim] - 1);
        stepLength = 1;
        QI::Log(verbose, "Fill dimension {} with gradient {}-{}", fillDim, startVal, stopVal);
    } else if (steps_arg) {
        fillType = FillTypes::Steps;
        std::istringstream stream(steps_arg.Get());
        stream >> fillDim;
        stream >> startVal;
        stream >> stopVal;
        int steps;
        stream >> steps;
        if (steps < 2) {
            QI::Fail("Must have more than 1 step");
        }
        stepLength = imgSize[fillDim] / steps;
        deltaVal   = (stopVal - startVal) / (steps - 1);
        QI::Log(verbose, "Fill dimension {} with {} steps {}-{}", fillDim, steps, startVal,
                stopVal);
    }
    QI::Log(verbose, "Step length = {} delta = {}", stepLength, deltaVal);
    imgRegion.SetIndex(imgIndex);
    imgRegion.SetSize(imgSize);
    newimg->SetRegions(imgRegion);
    newimg->SetSpacing(imgSpacing);
    newimg->SetOrigin(imgOrigin);
    newimg->Allocate();

    if (fillDim >= dim)
        QI::Fail("Fill dimension is larger than image dimension");
    itk::ImageLinearIteratorWithIndex<ImageType> it(newimg, imgRegion);
    it.SetDirection(fillDim);
    QI::Log(verbose, "Filling...");
    it.GoToBegin();
    while (!it.IsAtEnd()) {
        float val = startVal;
        while (!it.IsAtEndOfLine()) {
            it.Set(val);
            ++it;
            if (fillType == FillTypes::Gradient) {
                val += deltaVal;
            } else if ((it.GetIndex()[fillDim] % stepLength) == 0) {
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

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {

    QI::ParseArgs(parser, argc, argv, verbose);
    if (dims_arg.Get() == 3) {
        make_image<3>();
    } else if (dims_arg.Get() == 4) {
        make_image<4>();
    } else {
        std::cerr << "Unsupported dimension: " << dims_arg.Get() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
