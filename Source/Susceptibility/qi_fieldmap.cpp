/*
 *  qi_cestasym.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#include "Args.h"
#include "ImageIO.h"
#include "JSON.h"
#include "Util.h"

int fieldmap_main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Simple field-map via complex "
                                "division\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> input_path(parser, "INPUT", "Input multi-echo GRE file");
    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=hardware limit or $QUIT_THREADS)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(
        parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<double> delta_te(parser, "ΔTE", "Echo time difference (ms)", {"delta_te"});
    args::ValueFlag<double> B0(
        parser, "B0", "Field-strength in Tesla. Output will be in PPM", {"B0"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI::Log(verbose, "Opening file: {}", QI::CheckPos(input_path));
    auto input = QI::ReadImage<QI::VectorVolumeXF>(QI::CheckPos(input_path), verbose);
    QI::Log(verbose, "ΔTE = {} ms", QI::CheckValue(delta_te));

    auto fieldmap = QI::VolumeF::New();
    fieldmap->CopyInformation(input);
    fieldmap->SetRegions(input->GetBufferedRegion());
    fieldmap->Allocate(true);

    auto mt = itk::MultiThreaderBase::New();
    QI::Log(verbose, "Processing");
    const int N     = input->GetNumberOfComponentsPerPixel();
    auto      scale = 1e3 / (2. * M_PI * delta_te.Get()); // Convert to Hz
    if (B0) {
        const auto gamma = 42.57747892; // MHz per T to get PPM
        scale /= (gamma * B0.Get());
    }
    mt->SetNumberOfWorkUnits(threads.Get());
    mt->ParallelizeImageRegion<3>(
        input->GetBufferedRegion(),
        [&](const QI::VectorVolumeF::RegionType &region) {
            itk::ImageRegionConstIterator<QI::VectorVolumeXF> in_it(input, region);
            itk::ImageRegionIterator<QI::VolumeF>             f_it(fieldmap, region);
            for (in_it.GoToBegin(); !in_it.IsAtEnd(); ++in_it, ++f_it) {
                const Eigen::Map<const Eigen::ArrayXcf> data(in_it.Get().GetDataPointer(), 3);
                const auto phi = (data.tail(N - 1) / data.head(N - 1)).arg().mean();
                const auto f0  = phi * scale;
                f_it.Set(f0);
            }
        },
        nullptr);

    QI::WriteImage(fieldmap, outarg.Get() + "Fieldmap" + QI::OutExt(), verbose);
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
