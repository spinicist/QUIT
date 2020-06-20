/*
 *  qi_tvmask.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>
#include <iostream>

#include "Args.h"
#include "ImageIO.h"
#include "Util.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkMultiThreaderBase.h"

/*
 * Main
 */
int tvmask_main(args::Subparser &parser) {
    args::Positional<std::string> input_path(parser, "INPUT_FILE", "Input file");
    args::ValueFlag<int>          threads(parser,
                                 "THREADS",
                                 "Use N threads (default=hardware limit or $QUIT_THREADS)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string>  outarg(
        parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    args::ValueFlag<float> thresh(
        parser, "THRESH", "Threshold for TV mask (default 2)", {'t', "thresh"}, 2.0);

    parser.Parse();
    auto              input = QI::ReadImage<QI::VectorVolumeF>(QI::CheckPos(input_path), verbose);
    const std::string outPrefix = outarg ? outarg.Get() : QI::Basename(input_path.Get());
    const int         insize    = input->GetNumberOfComponentsPerPixel();

    auto output = QI::NewImageLike(input);

    QI::Log(verbose, "Processing");
    auto mt = itk::MultiThreaderBase::New();
    mt->ParallelizeImageRegion<3>(
        input->GetBufferedRegion(),
        [&](const QI::VectorVolumeF::RegionType &region) {
            itk::ImageRegionConstIterator<QI::VectorVolumeF> in_it(input, region);
            itk::ImageRegionIterator<QI::VolumeF>            out_it(output, region);

            for (in_it.GoToBegin(); !in_it.IsAtEnd(); ++in_it, ++out_it) {
                const Eigen::Map<const Eigen::ArrayXf> in_vec(in_it.Get().GetDataPointer(), insize);
                auto const tv    = (in_vec.tail(insize - 1) - in_vec.head(insize - 1)).abs().sum();
                auto const range = in_vec.maxCoeff() - in_vec.minCoeff();
                auto const ratio = tv / range;
                out_it.Set(ratio < thresh.Get());
            }
        },
        nullptr);
    const std::string fname = outPrefix + "tvmask" + QI::OutExt();
    QI::WriteImage(output, fname, verbose);
    return EXIT_SUCCESS;
}
