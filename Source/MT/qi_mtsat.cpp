/*
 *  qi_mtsat.cpp
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>

#include "Args.h"
#include "ImageIO.h"
#include "SequenceBase.h"
#include "Util.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

struct MTSatSequence : QI::SequenceBase {
    double TR, al_t1, al_pd, al_mt;

    QI_SEQUENCE_DECLARE(MTSat);
    Eigen::Index size() const override { return 1; };
};
void from_json(const json &j, MTSatSequence &s) {
    j.at("TR").get_to(s.TR);
    s.al_pd = j.at("FA_PD").get<double>() * M_PI / 180;
    s.al_t1 = j.at("FA_T1").get<double>() * M_PI / 180;
    s.al_mt = j.at("FA_MT").get<double>() * M_PI / 180;
}

void to_json(json &j, const MTSatSequence &s) {
    j = json{
        {"TR", s.TR},
        {"FA_PD", s.al_pd * 180 / M_PI},
        {"FA_T1", s.al_t1 * 180 / M_PI},
        {"FA_MT", s.al_mt * 180 / M_PI},
    };
}

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser(
        "Calculates MTsat from PDw, T1w, MTw data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> pdw_path(parser, "PDw", "Input PD-weighted file");
    args::Positional<std::string> t1w_path(parser, "T1w", "Input T1-weighted file");
    args::Positional<std::string> mtw_path(parser, "MTw", "Input MT-weighted file");

    args::HelpFlag       help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(
        parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    args::ValueFlag<std::string> mask_path(
        parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> json_file(
        parser, "JSON", "Read JSON from file instead of stdin", {"json"});
    args::ValueFlag<std::string> b1_path(parser, "B1", "Path to B1 map", {'b', "B1"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    auto pdw_img = QI::ReadImage(QI::CheckPos(pdw_path), verbose);
    auto t1w_img = QI::ReadImage(QI::CheckPos(t1w_path), verbose);
    auto mtw_img = QI::ReadImage(QI::CheckPos(mtw_path), verbose);
    QI::Log(verbose, "Reading sequence parameters");
    json          doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    MTSatSequence s(doc["MTSat"]);

    QI::VolumeF::Pointer b1_img = nullptr, mask_img = nullptr;
    if (b1_path) {
        b1_img = QI::ReadImage(b1_path.Get(), verbose);
    }
    if (mask_path) {
        mask_img = QI::ReadImage(mask_path.Get(), verbose);
    }

    QI::Info(verbose, "Allocating output memory");
    auto R1_img = QI::NewImageLike<QI::VolumeF>(pdw_img);
    auto A_img  = QI::NewImageLike<QI::VolumeF>(pdw_img);
    auto d_img  = QI::NewImageLike<QI::VolumeF>(pdw_img);

    QI::Info(verbose, "Processing");
    auto mt = itk::MultiThreaderBase::New();
    mt->SetNumberOfWorkUnits(threads.Get());
    mt->ParallelizeImageRegion<3>(
        R1_img->GetBufferedRegion(),
        [&](const QI::VolumeF::RegionType &region) {
            itk::ImageRegionConstIterator<QI::VolumeF> pdw_it(pdw_img, region),
                t1w_it(t1w_img, region), mtw_it(mtw_img, region);
            itk::ImageRegionIterator<QI::VolumeF> R1_it(R1_img, region), A_it(A_img, region),
                d_it(d_img, region);
            itk::ImageRegionConstIterator<QI::VolumeF> B1_it, mask_it;
            if (b1_img) {
                B1_it = itk::ImageRegionConstIterator<QI::VolumeF>(b1_img, region);
            }
            if (mask_img) {
                mask_it = itk::ImageRegionConstIterator<QI::VolumeF>(mask_img, region);
            }

            for (R1_it.GoToBegin(); !R1_it.IsAtEnd(); ++R1_it) {
                if (!mask_img || mask_it.Get()) {
                    auto S_pd = pdw_it.Get();
                    auto S_t1 = t1w_it.Get();
                    auto S_mt = mtw_it.Get();
                    auto B1   = b1_img ? B1_it.Get() : 1.0f;

                    auto R1 = (B1 * B1 / (2. * s.TR)) * (S_t1 * s.al_t1 - S_pd * s.al_pd) /
                              (S_pd / s.al_pd - S_t1 / s.al_t1);
                    auto A = (S_pd * S_t1 / B1) * (s.al_t1 / s.al_pd - s.al_pd / s.al_t1) /
                             (S_t1 * s.al_t1 - S_pd * s.al_pd);
                    auto d = (A * B1 * s.al_mt / S_mt - 1) * R1 * s.TR -
                             (B1 * B1 * s.al_mt * s.al_mt / 2);
                    R1_it.Set(R1);
                    A_it.Set(A);
                    d_it.Set(d * 100);
                } else {
                    R1_it.Set(0.0);
                    A_it.Set(0.0);
                    d_it.Set(0.0);
                }

                ++t1w_it;
                ++pdw_it;
                ++mtw_it;
                // ++R1_it;
                ++A_it;
                ++d_it;
                if (b1_img) {
                    ++B1_it;
                }
                if (mask_img) {
                    ++mask_it;
                }
            }
        },
        nullptr);
    QI::Info(verbose, "Finished");
    QI::WriteImage(R1_img, outarg.Get() + "MTSat_R1" + QI::OutExt(), verbose);
    QI::WriteImage(A_img, outarg.Get() + "MTSat_A" + QI::OutExt(), verbose);
    QI::WriteImage(d_img, outarg.Get() + "MTSat_delta" + QI::OutExt(), verbose);
    return EXIT_SUCCESS;
}
