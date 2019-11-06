/*
 *  qi_asl.cpp
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
#include "CASLSequence.h"
#include "ImageIO.h"
#include "Util.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

/*
 * Main
 */
int asl_main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates CBF from ASL data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "ASL_FILE", "Input ASL file");

    args::HelpFlag       help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string> json_file(
        parser, "JSON", "Read JSON from file instead of stdin", {"json"});
    args::ValueFlag<std::string> outarg(
        parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    args::ValueFlag<std::string> mask(
        parser, "MASK", "Only process voxels within the mask", {'m', "mask"});

    args::Flag average(parser, "AVERAGE", "Average the time-series", {'a', "average"});
    args::Flag slice_time(parser,
                          "SLICE TIME CORRECTION",
                          "Apply slice-time correction (number of post-label "
                          "delays must match number of slices)",
                          {'s', "slicetime"});

    args::ValueFlag<int> dummies(
        parser, "DUMMIES", "Discard this many dummy pairs from input", {'d', "dummies"}, 0);

    args::ValueFlag<std::string> T1_tissue_path(
        parser, "TISSUE T1", "Path to tissue T1 map (units are seconds)", {'t', "tissue"});
    args::ValueFlag<std::string> PD_path(parser, "PROTON DENSITY", "Path to PD image", {'p', "pd"});

    args::ValueFlag<double> T1_blood(parser,
                                     "BLOOD T1",
                                     "Value of blood T1 to use (seconds), default 1.65 for 3T",
                                     {'b', "blood"},
                                     1.65);
    args::ValueFlag<double> alpha(
        parser, "ALPHA", "Labelling efficiency, default 0.9", {'a', "alpha"}, 0.9);
    args::ValueFlag<double> lambda(parser,
                                   "LAMBDA",
                                   "Blood-brain partition co-efficent, default 0.9 mL/g",
                                   {'l', "lambda"},
                                   0.9);

    QI::ParseArgs(parser, argc, argv, verbose, threads);
    auto input = QI::ReadImage<QI::VectorVolumeF>(QI::CheckPos(input_path), verbose);
    QI::Log(verbose, "Reading sequence parameters");
    json                 doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    QI::CASLSequence     sequence(doc["CASL"]);
    QI::VolumeF::Pointer t1_img = nullptr, pd_img = nullptr, mask_img = nullptr;
    if (T1_tissue_path) {
        t1_img = QI::ReadImage(T1_tissue_path.Get(), verbose);
    }
    if (PD_path) {
        pd_img = QI::ReadImage(PD_path.Get(), verbose);
    }
    if (mask) {
        mask_img = QI::ReadImage(mask.Get(), verbose);
    }

    const auto n_slices = input->GetLargestPossibleRegion().GetSize()[2];
    if (slice_time && (n_slices != static_cast<size_t>(sequence.post_label_delay.rows()))) {
        QI::Fail("Number of post-label delays {} does not match number of slices {}",
                 sequence.post_label_delay.rows(),
                 n_slices);
    }

    const auto insize  = input->GetNumberOfComponentsPerPixel() - (dummies.Get() * 2);
    const auto outsize = average ? 1 : insize / 2;

    auto output = QI::VectorVolumeF::New();
    output->CopyInformation(input);
    output->SetRegions(input->GetBufferedRegion());
    output->SetNumberOfComponentsPerPixel(outsize);
    output->Allocate(true);
    auto const scales =
        6000 * lambda.Get() * exp(sequence.post_label_delay / T1_blood.Get()) /
        (2. * alpha.Get() * T1_blood.Get() * (1. - exp(-sequence.label_time / T1_blood.Get())));
    QI::Info(verbose, "Processing");
    auto mt = itk::MultiThreaderBase::New();
    mt->SetNumberOfWorkUnits(threads.Get());
    mt->ParallelizeImageRegion<3>(
        input->GetBufferedRegion(),
        [&](const QI::VectorVolumeF::RegionType &region) {
            itk::ImageRegionConstIterator<QI::VectorVolumeF>     in_it(input, region);
            itk::ImageRegionIteratorWithIndex<QI::VectorVolumeF> out_it(output, region);
            itk::ImageRegionConstIterator<QI::VolumeF>           t1_it, pd_it, mask_it;
            if (t1_img) {
                t1_it = itk::ImageRegionConstIterator<QI::VolumeF>(t1_img, region);
                t1_it.GoToBegin();
            }
            if (pd_img) {
                pd_it = itk::ImageRegionConstIterator<QI::VolumeF>(pd_img, region);
                pd_it.GoToBegin();
            }
            if (mask_img) {
                mask_it = itk::ImageRegionConstIterator<QI::VolumeF>(mask_img, region);
                mask_it.GoToBegin();
            }
            out_it.GoToBegin();
            for (in_it.GoToBegin(); !in_it.IsAtEnd(); ++in_it) {
                itk::VariableLengthVector<float> CBF_vector(outsize);
                if (!mask_img || mask_it.Get()) {
                    const Eigen::Map<const Eigen::ArrayXXf> raw(
                        in_it.Get().GetDataPointer() + (dummies.Get() * 2), 2, insize / 2);
                    const auto &label        = raw.row(0);
                    const auto &control      = raw.row(1);
                    const auto &difference   = control - label; // Negative contrast
                    const auto PD_correction = t1_img ? 1.0 - exp(-sequence.TR / t1_it.Get()) : 1.0;
                    const auto PD    = (pd_img ? pd_it.Get() : control.mean()) / PD_correction;
                    const auto scale = sequence.post_label_delay.rows() > 1 ?
                                           scales[out_it.GetIndex()[2]] :
                                           scales[0];
                    const auto CBF = scale * difference / PD;
                    if (average) {
                        CBF_vector[0] = CBF.mean();
                    } else {
                        for (size_t i = 0; i < outsize; i++) {
                            CBF_vector[i] = CBF[i];
                        }
                    }
                } else {
                    CBF_vector.Fill(0.0);
                }
                out_it.Set(CBF_vector);
                ++out_it;
                if (t1_img) {
                    ++t1_it;
                }
                if (pd_img) {
                    ++pd_it;
                }
                if (mask_img) {
                    ++mask_it;
                }
            }
        },
        nullptr);
    QI::Info(verbose, "Finished");
    QI::WriteImage(output, outarg.Get() + "CASL_CBF" + QI::OutExt(), verbose);
    return EXIT_SUCCESS;
}
