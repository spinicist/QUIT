/*
 *  qi_b1_papp.cpp
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <string>

#include "Args.h"
#include "ImageIO.h"
#include "ImageTypes.h"
#include "Util.h"

#include "itkDivideImageFilter.h"
#include "itkExtractImageFilter.h"

int main(int argc, char **argv) {
    args::ArgumentParser parser("Calculates B1- (receive bias) maps. Input file should have two "
                                "volumes - body coil then head coil\n"
                                "http://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "INPUT", "Input file");
    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string> out_prefix(
        parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    auto inFile    = QI::ReadImage<QI::SeriesF>(QI::CheckPos(input_path), verbose);
    auto body_coil = itk::ExtractImageFilter<QI::SeriesF, QI::VolumeF>::New();
    auto head_coil = itk::ExtractImageFilter<QI::SeriesF, QI::VolumeF>::New();
    auto region    = inFile->GetLargestPossibleRegion();

    region.GetModifiableSize()[3]  = 0;
    region.GetModifiableIndex()[3] = 0;
    body_coil->SetExtractionRegion(region);
    body_coil->SetInput(inFile);
    body_coil->SetDirectionCollapseToSubmatrix();
    region.GetModifiableIndex()[3] = 1;
    head_coil->SetExtractionRegion(region);
    head_coil->SetInput(inFile);
    head_coil->SetDirectionCollapseToSubmatrix();

    auto ratio = itk::DivideImageFilter<QI::VolumeF, QI::VolumeF, QI::VolumeF>::New();
    ratio->SetInput(0, head_coil->GetOutput());
    ratio->SetInput(1, body_coil->GetOutput());
    QI::WriteImage(ratio->GetOutput(), out_prefix.Get() + "B1minus" + QI::OutExt(), verbose);
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
