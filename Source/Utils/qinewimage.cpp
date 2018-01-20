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

#include <iostream>
#include <random>

#include "itkImage.h"
#include "itkImageSliceIteratorWithIndex.h"

#include "Util.h"
#include "IO.h"
#include "Args.h"

/*
 * Define arguments globally because I am lazy
 */
args::ArgumentParser parser("This is a tool to create 3D Nifti files, either blank headers with orientation "
                            "information, e.g. for registration, or files filled with simple patterns of "
                            "data e.g. solid values, gradients, or blocks. The default is to create a 3D file "
                            "filled with zeros. Choose from the options below to create something else."
                            "\nhttp://github.com/spinicist/QUIT");

args::Positional<std::string> fName(parser, "OUTPUT", "Output filename");

args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
args::ValueFlag<int> dims_arg(parser, "DIMS", "Image dimension, default 3", {'d', "dims"}, 3);
args::ValueFlag<std::string> size_arg(parser, "SIZE", "Image size", {'s', "size"});
args::ValueFlag<std::string> spacing_arg(parser, "SPACING", "Voxel spacing", {'p', "spacing"});
args::ValueFlag<std::string> origin_arg(parser, "ORIGIN", "Image origin", {'o', "origin"});
args::ValueFlag<std::string> fill_arg(parser, "FILL", "Fill with value", {'f', "fill"});
args::ValueFlag<std::string> grad_arg(parser, "GRADIENT", "Fill with gradient (dim, low, high)", {'g', "grad"});
args::ValueFlag<std::string> steps_arg(parser, "STEPS", "Fill with discrete steps (dim, low, high, steps)", {'t', "step"});

template<int dim>
void make_image() {
    typedef itk::Image<float, dim> ImageType;
    using RegionType = typename ImageType::RegionType;
    using IndexType = typename ImageType::IndexType;
    using SizeType = typename ImageType::SizeType;
    using SpacingType = typename ImageType::SpacingType;
    using PointType = typename ImageType::PointType;
    auto newimg = ImageType::New();
    RegionType  imgRegion;
    IndexType   imgIndex;   imgIndex.Fill(0);
    SizeType    imgSize;    imgSize.Fill(1);
    SpacingType imgSpacing; imgSpacing.Fill(1.0);
    PointType   imgOrigin;  imgOrigin.Fill(0.);
    
    if (size_arg)    QI::ArrayArg<SizeType, dim>(size_arg.Get(), imgSize);
    if (spacing_arg) QI::ArrayArg<SpacingType, dim>(spacing_arg.Get(), imgSpacing);
    if (origin_arg)  QI::ArrayArg<PointType, dim>(origin_arg.Get(), imgOrigin);

    if (verbose) {
        std::cout << "Size:    " << imgSize << std::endl;
        std::cout << "Spacing: " << imgSpacing << std::endl;
        std::cout << "Origin:  " << imgOrigin << std::endl;
    }
    enum class FillTypes { Fill, Gradient, Steps };
    FillTypes fillType = FillTypes::Fill;
    float startVal = 0, deltaVal = 0, stopVal = 0;
    int stepLength = 1, fillDim = 0;
    if (fill_arg) {
        fillType = FillTypes::Fill;
        startVal = std::stod(fill_arg.Get());
        if (verbose) std::cout << "Fill with constant value: " << startVal << std::endl;
    } else if (grad_arg) {
        fillType = FillTypes::Gradient;
        std::istringstream stream(grad_arg.Get());
        stream >> fillDim;
        stream >> startVal;
        stream >> stopVal;
        deltaVal = (stopVal - startVal) / (imgSize[fillDim] - 1);
        stepLength = 1;
        if (verbose) std::cout << "Fill with gradient " << startVal << "-" << stopVal << " on dim " << fillDim << std::endl;
    } else if (steps_arg) {
        fillType = FillTypes::Steps;
        std::istringstream stream(steps_arg.Get());
        stream >> fillDim;
        stream >> startVal;
        stream >> stopVal;
        int steps;
        stream >> steps;
        if (steps < 2) {
            QI_FAIL("Must have more than 1 step");
        }
        stepLength = imgSize[fillDim] / steps;
        deltaVal = (stopVal - startVal) / (steps - 1);
        if (verbose) std::cout << "Fill with " << steps << " steps " << startVal << "-" << stopVal << " on dim " << fillDim << std::endl;
    }
    if (verbose) std::cout << "Step length is " << stepLength << " delta value is " << deltaVal << std::endl;
    imgRegion.SetIndex(imgIndex);
    imgRegion.SetSize(imgSize);
    newimg->SetRegions(imgRegion);
    newimg->SetSpacing(imgSpacing);
    newimg->SetOrigin(imgOrigin);
    newimg->Allocate();

    if (fillDim >= dim) QI_FAIL("Fill dimension is larger than image dimension");
    itk::ImageSliceIteratorWithIndex<ImageType> it(newimg, imgRegion);
    it.SetFirstDirection((fillDim + 1) % dim);
    it.SetSecondDirection((fillDim + 2) % dim);
    float val = startVal;
    if (verbose) std::cout << "Filling..." << std::endl;
    it.GoToBegin();
    while (!it.IsAtEnd()) {
        while (!it.IsAtEndOfSlice()) {
            while (!it.IsAtEndOfLine()) {
                it.Set(val);
                ++it;
            }
            it.NextLine();
        }
        it.NextSlice();
        if (fillType == FillTypes::Gradient) {
            val += deltaVal;
        } else if ((it.GetIndex()[fillDim] % stepLength) == 0) {
            val += deltaVal;
        }
    }
    if (verbose) std::cout << "Writing file to: " << QI::CheckPos(fName) << std::endl;
    QI::WriteImage(newimg, QI::CheckPos(fName));
}

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {

    QI::ParseArgs(parser, argc, argv);
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


