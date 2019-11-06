/*
 *  qi_select.cpp
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

#include "itkExtractImageFilter.h"
#include "itkTileImageFilter.h"

int select_main(int argc, char **argv) {
    args::ArgumentParser parser(
        "Extracts a set of volumes from a time-series and saves as a new series\n"
        "http://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "INPUT", "Input file");
    args::Positional<std::string> output_path(parser, "OUTPUT", "Output file");
    args::Positional<std::string> volume_list(parser, "VOLUMES", "Comma separated list of volumes");
    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    auto in_file        = QI::ReadImage<QI::SeriesF>(QI::CheckPos(input_path), verbose);
    auto volume_indices = QI::IntsFromString(volume_list.Get());

    auto tiler  = itk::TileImageFilter<QI::VolumeF, QI::SeriesF>::New();
    auto region = in_file->GetLargestPossibleRegion();

    itk::FixedArray<unsigned int, 4> layout;
    layout[0] = layout[1] = layout[2] = 1;
    layout[3]                         = volume_indices.size();
    tiler->SetLayout(layout);

    int const max_index           = region.GetSize()[3] - 1;
    region.GetModifiableSize()[3] = 0;

    using Extractor = itk::ExtractImageFilter<QI::SeriesF, QI::VolumeF>;
    std::vector<Extractor::Pointer> volumes(volume_indices.size());

    for (auto ii = 0; ii < static_cast<int>(volume_indices.size()); ii++) {
        auto const &volume_index = volume_indices[ii];
        if (volume_index > max_index) {
            QI::Fail("Invalid volume index {}, max is {}", volume_index, max_index);
        }
        QI::Info(verbose, "Out volume {} = in volume {}", ii, volume_index);
        region.GetModifiableIndex()[3] = volume_index;
        auto volume                    = itk::ExtractImageFilter<QI::SeriesF, QI::VolumeF>::New();
        volume->SetExtractionRegion(region);
        volume->SetInput(in_file);
        volume->SetDirectionCollapseToSubmatrix();
        tiler->SetInput(ii, volume->GetOutput());
        volumes.emplace_back(volume);
    }
    tiler->Update();
    QI::WriteImage(tiler->GetOutput(), output_path.Get(), verbose);
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
