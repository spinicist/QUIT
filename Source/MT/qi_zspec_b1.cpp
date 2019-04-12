/*
 *  qi_cestasym.cpp
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>
#include <algorithm>
#include <iterator>

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

// #define QI_DEBUG_BUILD 1

#include "Args.h"
#include "ImageIO.h"
#include "JSON.h"
#include "Macro.h"
#include "Spline.h"
#include "Util.h"

int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser(
        "B1-correction for Z-spectrums using linear "
        "interpolation as in Windschuh et al\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string>     b1_path(parser, "B1", "Path to B1-map");
    args::PositionalList<std::string> input_paths(
        parser, "INPUTS", "Input Z-spectra files (1 file per B1 level)");
    args::HelpFlag       help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string> out_suffix(
        parser, "OUTPUT", "Change ouput filenames (default is input_b1)", {'o', "out"}, "_b1");
    args::ValueFlag<std::string> json_file(
        parser, "JSON", "Read JSON from file instead of stdin", {"json"});
    args::ValueFlag<std::string> mask(
        parser, "MASK", "Only process voxels within the mask", {'m', "mask"});

    QI::ParseArgs(parser, argc, argv, verbose, threads);

    const auto b1_image = QI::ReadImage(QI::CheckPos(b1_path), verbose);
    QI::CheckList(input_paths);
    std::vector<QI::VectorVolumeF::Pointer> inputs, outputs;
    for (auto const &ip : input_paths.Get()) {
        inputs.emplace_back(QI::ReadImage<QI::VectorVolumeF>(ip, verbose));
        outputs.emplace_back(QI::VectorVolumeF::New());
        outputs.back()->CopyInformation(b1_image);
        outputs.back()->SetRegions(b1_image->GetBufferedRegion());
        outputs.back()->SetNumberOfComponentsPerPixel(
            inputs.back()->GetNumberOfComponentsPerPixel());
        outputs.back()->Allocate(true);
    }
    const Eigen::Index Nz = inputs.front()->GetNumberOfComponentsPerPixel();

    for (auto const &i : inputs) {
        if (i->GetNumberOfComponentsPerPixel() != Nz) {
            QI::Fail("All input Z-spectra must be the same length");
        }
    }

    json doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);

    auto b1_rms = QI::ArrayFromJSON(doc, "b1_rms");
    if (b1_rms.rows() != static_cast<Eigen::Index>(inputs.size())) {
        QI::Fail("The number of B1 RMS entries must match the number of inputs");
    }
    QI::Log(verbose, "B1 RMS Powers: {}", b1_rms.transpose());

    const QI::VolumeF::Pointer mask_image = mask ? QI::ReadImage(mask.Get(), verbose) : nullptr;

    auto mt = itk::MultiThreaderBase::New();
    QI::Log(verbose, "Processing");
    mt->SetNumberOfWorkUnits(threads.Get());
    mt->ParallelizeImageRegion<3>(
        inputs.front()->GetBufferedRegion(),
        [&](const QI::VectorVolumeF::RegionType &region) {
            auto b1_it = itk::ImageRegionConstIterator<QI::VolumeF>(b1_image, region);

            std::vector<itk::ImageRegionConstIterator<QI::VectorVolumeF>> in_its;
            std::vector<itk::ImageRegionIterator<QI::VectorVolumeF>>      out_its;
            for (auto const &i : inputs) {
                in_its.emplace_back(i, region);
            }
            for (auto const &o : outputs) {
                out_its.emplace_back(o, region);
            }

            itk::ImageRegionConstIterator<QI::VolumeF> mask_it;
            if (mask_image)
                mask_it = itk::ImageRegionConstIterator<QI::VolumeF>(mask_image, region);

            for (b1_it.GoToBegin(); !b1_it.IsAtEnd(); ++b1_it) {
                const double B1 = b1_it.Get();
                if (!mask_image || mask_it.Get()) {
                    std::vector<const Eigen::Map<const Eigen::ArrayXf>> in_z;
                    std::vector<itk::VariableLengthVector<float>>       out_z;
                    for (const auto &ii : in_its) {
                        in_z.emplace_back(ii.Get().GetDataPointer(), Nz);
                        out_z.emplace_back(Nz);
                    }

                    for (Eigen::Index iz = 0; iz < Nz; iz++) {
                        Eigen::VectorXd Z(b1_rms.rows());
                        for (Eigen::Index ib = 0; ib < b1_rms.rows(); ib++) {
                            Z[ib] = in_z[ib][iz];
                        }
                        Eigen::VectorXd slope = (Z.transpose() * Z)
                                                    .partialPivLu()
                                                    .solve(Z.transpose() * B1 * b1_rms.matrix());
                        Eigen::VectorXd Zcorrected = slope[0] * b1_rms;
                        for (Eigen::Index ib = 0; ib < b1_rms.rows(); ib++) {
                            out_z[ib][iz] = Zcorrected[ib];
                        }
                    }

                    for (size_t i = 0; i < out_its.size(); i++) {
                        out_its[i].Set(out_z[i]);
                    }
                } else {
                    itk::VariableLengthVector<float> zero(Nz);
                    zero.Fill(0.0);
                    for (auto &oi : out_its) {
                        oi.Set(zero);
                    }
                }

                for (auto &ii : in_its) {
                    ++ii;
                }
                for (auto &oi : out_its) {
                    ++oi;
                }
                if (mask_image)
                    ++mask_it;
            }
        },
        nullptr);

    for (size_t io = 0; io < outputs.size(); io++) {
        std::string outname =
            QI::StripExt(QI::Basename(input_paths.Get()[io])) + out_suffix.Get() + QI::OutExt();
        QI::WriteImage(outputs[io], outname, verbose);
    }
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
