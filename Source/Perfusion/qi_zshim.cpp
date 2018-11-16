/*
 *  qi_zshim.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Core>

#include "Util.h"
#include "ImageIO.h"
#include "Args.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Combines Z-Shimmed data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "ZSHIM_FILE", "Input Z-Shimmed file");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Output progress messages", {'v', "verbose"});
    args::Flag     debug(parser, "DEBUG", "Output debug messages", {'d', "debug"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    args::ValueFlag<int> zshims(parser, "ZSHIMS", "Number of Z-Shims (default 8)", {'z', "zshims"}, 8);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    auto input = QI::ReadVectorImage(QI::CheckPos(input_path), verbose);
    const std::string outPrefix = outarg ? outarg.Get() : QI::Basename(input_path.Get());
    const auto insize  = input->GetNumberOfComponentsPerPixel();
    const auto outsize = insize / zshims.Get();

    auto output = QI::VectorVolumeF::New();
    output->CopyInformation(input);
    output->SetRegions(input->GetBufferedRegion());
    output->SetNumberOfComponentsPerPixel(outsize);
    output->Allocate(true);

    QI_LOG(verbose, "Processing");
    auto mt = itk::MultiThreaderBase::New();
    mt->ParallelizeImageRegion<3>(input->GetBufferedRegion(),
        [&](const QI::VectorVolumeF::RegionType &region) {
            itk::ImageRegionConstIterator<QI::VectorVolumeF> in_it(input, region);
            itk::ImageRegionIterator<QI::VectorVolumeF> out_it(output, region);

            for (in_it.GoToBegin(); !in_it.IsAtEnd(); ++in_it, ++out_it) {
                const Eigen::Map<const Eigen::MatrixXf> input(in_it.Get().GetDataPointer(), zshims.Get(), outsize);
                // std::cout << "***\n" << input << "\n" << std::endl;
                const Eigen::VectorXf shimmed = input.colwise().norm();
                // std::cout << shimmed.transpose() << std::endl;
                itk::VariableLengthVector<float> itk_shimmed(shimmed.data(), outsize);
                out_it.Set(itk_shimmed);
            }
        },
        nullptr);
    const std::string fname = outPrefix + "_zshim" + QI::OutExt();
    QI_LOG(verbose, "Writing file: " << fname);
    QI::WriteVectorImage(output, fname);
    return EXIT_SUCCESS;
}
