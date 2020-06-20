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

#include <list>
#include <memory>

#include "Args.h"
#include "ImageIO.h"
#include "ImageTypes.h"
#include "PathUnwrapFilter.h"
#include "ReliabilityFilter.h"
#include "Util.h"
#include "itkExtractImageFilter.h"
#include "itkImageToImageFilter.h"
#include "itkTileImageFilter.h"

/*
 * Main
 */
int unwrap_path_main(args::Subparser &parser) {
    args::Positional<std::string> input_path(parser, "PHASE", "Wrapped phase image");
    args::ValueFlag<int>          threads(parser,
                                 "THREADS",
                                 "Use N threads (default=hardware limit or $QUIT_THREADS)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string>  outarg(
        parser, "OUTPUT PREFIX", "Change output prefix (default input filename)", {'o', "out"});
    args::ValueFlag<std::string> maskarg(
        parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    parser.Parse();

    auto inFile = QI::ReadImage<QI::SeriesF>(QI::CheckPos(input_path), verbose);

    typedef itk::ExtractImageFilter<QI::SeriesF, QI::VolumeF> TExtract;
    typedef itk::TileImageFilter<QI::VolumeF, QI::SeriesF>    TTile;

    auto         region           = inFile->GetLargestPossibleRegion();
    const size_t nvols            = region.GetSize()[3]; // Save for the loop
    region.GetModifiableSize()[3] = 0;
    itk::FixedArray<unsigned int, 4> layout;
    layout[0] = layout[1] = layout[2] = 1;
    layout[3]                         = nvols;

    auto extract = TExtract::New();
    auto tile    = TTile::New();
    extract->SetInput(inFile);
    extract->SetDirectionCollapseToSubmatrix();
    tile->SetLayout(layout);

    auto reliabilityFilter = itk::PhaseReliabilityFilter::New();
    auto unwrapFilter      = itk::UnwrapPathPhaseFilter::New();
    for (size_t i = 0; i < nvols; i++) {
        region.GetModifiableIndex()[3] = i;
        QI::Log(verbose, "Processing volume {}", i);
        extract->SetExtractionRegion(region);
        extract->Update();
        QI::Log(verbose, "Calculating reliabilty");
        reliabilityFilter->SetInput(extract->GetOutput());
        reliabilityFilter->Update();
        QI::Log(verbose, "Unwrapping phase");
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
    inFile   = tile->GetOutput();
    inFile->SetDirection(dir);
    inFile->SetSpacing(spc);
    inFile->DisconnectPipeline();

    std::string outname =
        (outarg ? outarg.Get() : (QI::StripExt(input_path.Get())) + "_unwrapped" + QI::OutExt());
    QI::WriteImage(inFile, outname, verbose);

    return EXIT_SUCCESS;
}
