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
#include "Spline.h"
#include "Util.h"

int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Interpolates Z-spectrums using Cubic "
                                "Splines\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> input_path(parser, "INPUT", "Input Z-spectrum file");
    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(
        parser, "OUTPUT", "Change ouput filename (default is input_interp)", {'o', "out"});
    args::ValueFlag<std::string> json_file(
        parser, "JSON", "Read JSON from file instead of stdin", {"json"});
    args::ValueFlag<std::string> f0_arg(parser,
                                        "OFF RESONANCE",
                                        "Specify off-resonance frequency (units must match input)",
                                        {'f', "f0"});
    args::ValueFlag<int>         order(
        parser, "ORDER", "Spline order for interpolation (default 3)", {'O', "order"}, 3);
    args::Flag asym(parser, "ASYMMETRY", "Output asymmetry (M+ - M-)", {'a', "asym"});
    args::ValueFlag<std::string> ref_arg(
        parser,
        "REFERENCE",
        "Divide output by reference and multiply by 100 (output %)",
        {'r', "ref"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI::Log(verbose, "Opening file: {}", QI::CheckPos(input_path));
    auto input = QI::ReadImage<QI::VectorVolumeF>(QI::CheckPos(input_path), verbose);

    QI::Log(verbose, "Reading input frequencies");
    rapidjson::Document json = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    auto                in_freqs  = QI::ArrayFromJSON(json, "input_freqs");
    auto                out_freqs = QI::ArrayFromJSON(json, "output_freqs");
    QI::Log(verbose, "Input frequencies: {}", in_freqs.transpose());
    QI::Log(verbose, "Output frequencies:{}", out_freqs.transpose());
    if (asym)
        QI::Log(verbose, "Asymmetry output selected");
    if (order)
        QI::Log(verbose, "Spline order: {}", order.Get());

    const QI::VolumeF::Pointer f0_image = f0_arg ? QI::ReadImage(f0_arg.Get(), verbose) : nullptr;
    const QI::VolumeF::Pointer ref_image =
        ref_arg ? QI::ReadImage(ref_arg.Get(), verbose) : nullptr;
    QI::VectorVolumeF::Pointer output = QI::VectorVolumeF::New();
    output->CopyInformation(input);
    output->SetRegions(input->GetBufferedRegion());
    output->SetNumberOfComponentsPerPixel(out_freqs.rows());
    output->Allocate(true);

    auto mt = itk::MultiThreaderBase::New();
    QI::Log(verbose, "Processing");
    mt->SetNumberOfWorkUnits(threads.Get());
    mt->ParallelizeImageRegion<3>(
        input->GetBufferedRegion(),
        [&](const QI::VectorVolumeF::RegionType &region) {
            itk::ImageRegionConstIterator<QI::VectorVolumeF> in_it(input, region);
            itk::ImageRegionConstIterator<QI::VolumeF>       f0_it, ref_it;
            if (f0_arg)
                f0_it = itk::ImageRegionConstIterator<QI::VolumeF>(f0_image, region);
            if (ref_arg)
                ref_it = itk::ImageRegionConstIterator<QI::VolumeF>(ref_image, region);
            itk::ImageRegionIterator<QI::VectorVolumeF> out_it(output, region);
            for (in_it.GoToBegin(); !in_it.IsAtEnd(); ++in_it, ++out_it) {
                const Eigen::Map<const Eigen::ArrayXf> zdata(in_it.Get().GetDataPointer(),
                                                             in_freqs.rows());
                QI::SplineInterpolator           zspec(in_freqs, zdata.cast<double>(), order.Get());
                itk::VariableLengthVector<float> interped(out_freqs.rows());
                float                            f0 = f0_image ? f0_it.Get() : 0.0;
                for (int f = 0; f < out_freqs.rows(); f++) {
                    if (asym) {
                        const double p_frq = f0 + out_freqs[f];
                        const double m_frq = f0 - out_freqs[f];
                        interped[f]        = zspec(m_frq) - zspec(p_frq);
                    } else {
                        const double frq = f0 + out_freqs[f];
                        interped[f]      = zspec(frq);
                    }
                }
                if (ref_arg) {
                    interped *= (100. / ref_it.Get());
                    ++ref_it;
                }
                out_it.Set(interped);
                if (f0_image)
                    ++f0_it;
            }
        },
        nullptr);
    std::string outname =
        outarg ? outarg.Get() : QI::StripExt(input_path.Get()) + "_interp" + QI::OutExt();
    QI::WriteImage(output, outname, verbose);
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
