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

#include <Eigen/Core>

#include "Args.h"
#include "ImageIO.h"
#include "Util.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

/*
 * Main
 */
int zshim_main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Combines Z-Shimmed data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "ZSHIM_FILE", "Input Z-Shimmed file");
    args::HelpFlag                help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Output progress messages", {'v', "verbose"});
    args::Flag           debug(parser, "DEBUG", "Output debug messages", {'d', "debug"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=hardware limit or $QUIT_THREADS)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(
        parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    args::ValueFlag<int> zshims(
        parser, "ZSHIMS", "Number of Z-Shims (default 8)", {'z', "zshims"}, 8);
    args::ValueFlag<int> yshims(
        parser, "YSHIMS", "Number of Y-Shims (default 1)", {'y', "yshims"}, 1);
    args::ValueFlag<int>         zdrop(parser, "ZDROP", "Z-Shims to drop, default 0", {"zdrop"}, 0);
    args::ValueFlag<int>         ydrop(parser, "YDROP", "Y-Shims to drop, default 0", {"ydrop"}, 0);
    args::ValueFlag<std::string> noise(
        parser, "NOISE REGION", "Subtract noise measured in region", {'n', "noiseregion"});

    QI::ParseArgs(parser, argc, argv, verbose, threads);
    auto              input = QI::ReadImage<QI::VectorVolumeF>(QI::CheckPos(input_path), verbose);
    const std::string outPrefix = outarg ? outarg.Get() : QI::Basename(input_path.Get());
    const int         insize    = input->GetNumberOfComponentsPerPixel();
    const auto        gridsize  = zshims.Get() * yshims.Get();
    const auto        outsize   = insize / gridsize;

    auto output = QI::VectorVolumeF::New();
    output->CopyInformation(input);
    output->SetRegions(input->GetBufferedRegion());
    output->SetNumberOfComponentsPerPixel(outsize);
    output->Allocate(true);

    double noise_mean = 0., noise_sqr_mean = 0., noise_sigma = 0.;
    if (noise) {
        QI::Log(verbose, "Calculating noise statistics");
        auto mt           = itk::MultiThreaderBase::New();
        auto noise_region = QI::RegionFromString<QI::VectorVolumeF::RegionType>(noise.Get());
        itk::ImageRegionConstIterator<QI::VectorVolumeF> noise_it(input, noise_region);
        long                                             samples = 0;
        for (noise_it.GoToBegin(); !noise_it.IsAtEnd(); ++noise_it) {
            Eigen::Map<Eigen::VectorXf const> const noise_vec(noise_it.Get().GetDataPointer(),
                                                              insize);
            for (Eigen::Index ii = 0; ii < insize; ii++) {
                auto const val = noise_vec[ii];
                noise_mean += val;
                noise_sqr_mean += val * val;
                samples++;
            }
        }
        noise_mean     = noise_mean / samples;
        noise_sqr_mean = noise_sqr_mean / samples;
        noise_sigma    = sqrt(noise_sqr_mean - noise_mean * noise_mean);
        QI::Log(verbose,
                "Noise samples {} mean {:.2f} sd {:.2f} ratio {:.2f} squared mean {:.2f}",
                samples,
                noise_mean,
                noise_sigma,
                noise_mean / noise_sigma,
                noise_sqr_mean);
    }

    QI::Log(verbose, "Processing");
    auto mt = itk::MultiThreaderBase::New();
    mt->SetNumberOfWorkUnits(threads.Get());
    mt->ParallelizeImageRegion<3>(
        input->GetBufferedRegion(),
        [&](const QI::VectorVolumeF::RegionType &region) {
            itk::ImageRegionConstIterator<QI::VectorVolumeF> in_it(input, region);
            itk::ImageRegionIterator<QI::VectorVolumeF>      out_it(output, region);

            for (in_it.GoToBegin(); !in_it.IsAtEnd(); ++in_it, ++out_it) {
                itk::VariableLengthVector<float> itk_shimmed(outsize);
                for (auto out = 0; out < outsize; out++) {
                    const Eigen::Map<const Eigen::ArrayXf> grid(
                        in_it.Get().GetDataPointer() + gridsize * out, zshims.Get(), yshims.Get());
                    double val = (grid.block(zdrop.Get(),
                                             ydrop.Get(),
                                             zshims.Get() - 2 * zdrop.Get(),
                                             yshims.Get() - 2 * ydrop.Get())
                                      .square() -
                                  noise_sqr_mean)
                                     .sum();
                    itk_shimmed[out] = sqrt(std::max(val, 0.));
                }
                out_it.Set(itk_shimmed);
            }
        },
        nullptr);
    const std::string fname = outPrefix + "_zshim" + QI::OutExt();
    QI::WriteImage(output, fname, verbose);
    return EXIT_SUCCESS;
}
