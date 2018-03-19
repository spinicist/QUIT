/*
 *  qi_unwrap_path.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  This is an implementation of the algorithm found in:
 *  Abdul-Rahman et al, Fast and robust three-dimensional best path phase unwrapping algorithm,
 *  http://ao.osa.org/abstract.cfm?URI=ao-46-26-6623
 *  Abdul-Rahman et al, Robust three-dimensional best-path phase-unwrapping algorithm that 
 *  avoids singularity loops, http://ao.osa.org/abstract.cfm?URI=ao-48-23-4582
 */

#include <iostream>
#include <memory>
#include <list>
#include <list>

#include "itkImageToImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkTileImageFilter.h"
#include "ImageTypes.h"
#include "Util.h"
#include "ImageIO.h"
#include "Args.h"
#include "ReliabilityFilter.h"
#include "PathUnwrapFilter.h"

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Path-based phase unwrapping\n"
                                "See Abdul-Rahman et al. Fast and Robust three-dimensional best path phase unwrapping algorithm\n"
                                "http://ao.osa.org/abstract.cfm?URI=ao-46-26-6623\n"
                                "http://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "PHASE", "Wrapped phase image");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPUT PREFIX", "Change output prefix (default input filename)", {'o', "out"});
    args::ValueFlag<std::string> maskarg(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    QI::ParseArgs(parser, argc, argv);
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(threads.Get());

    if (verbose) std::cout << "Reading phase file: " << QI::CheckPos(input_path) << std::endl;
    auto inFile = QI::ReadImage<QI::SeriesF>(QI::CheckPos(input_path));

    typedef itk::ExtractImageFilter<QI::SeriesF, QI::VolumeF> TExtract;
    typedef itk::TileImageFilter<QI::VolumeF, QI::SeriesF>    TTile;
    
    auto region = inFile->GetLargestPossibleRegion();
    const size_t nvols = region.GetSize()[3]; // Save for the loop
    region.GetModifiableSize()[3] = 0;
    itk::FixedArray<unsigned int, 4> layout;
    layout[0] = layout[1] = layout[2] = 1;
    layout[3] = nvols;

    auto extract = TExtract::New();
    auto tile   = TTile::New();
    extract->SetInput(inFile);
    extract->SetDirectionCollapseToSubmatrix();
    tile->SetLayout(layout);

    auto reliabilityFilter = itk::PhaseReliabilityFilter::New();
    auto unwrapFilter = itk::UnwrapPathPhaseFilter::New();
    for (int i = 0; i < nvols; i++) {
        region.GetModifiableIndex()[3] = i;
        if (verbose) std::cout << "Processing volume " << i << std::endl;
        extract->SetExtractionRegion(region);
        extract->Update();
        if (verbose) std::cout << "Calculating reliabilty" << std::endl;
        reliabilityFilter->SetInput(extract->GetOutput());
        reliabilityFilter->Update();
        if (verbose) std::cout << "Unwrapping phase" << std::endl;
        unwrapFilter->SetInput(extract->GetOutput());
        unwrapFilter->SetReliability(reliabilityFilter->GetOutput());
        unwrapFilter->Update();
        tile->SetInput(i, unwrapFilter->GetOutput());
        unwrapFilter->GetOutput()->DisconnectPipeline();
    }
    tile->Update();
    // Make sure output orientation info is correct
    auto dir = inFile->GetDirection();
    auto spc = inFile->GetSpacing();
    inFile = tile->GetOutput();
    inFile->SetDirection(dir);
    inFile->SetSpacing(spc);
    inFile->DisconnectPipeline();

    std::string outname = (outarg ? outarg.Get() : (QI::StripExt(input_path.Get())) + "_unwrapped" + QI::OutExt());
    if (verbose) std::cout << "Writing output: " << outname << std::endl;
    QI::WriteImage(inFile, outname);

    return EXIT_SUCCESS;
}
