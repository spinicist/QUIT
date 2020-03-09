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

struct MTContrast {
    std::string               name;
    std::vector<Eigen::Index> add_indices, sub_indices, ref_indices;
    float                     scale   = 1.f;
    bool                      reverse = false, inverse = false;
};

void from_json(const json &j, MTContrast &s) {
    j.at("name").get_to(s.name);
    j.at("add").get_to(s.add_indices);
    j.at("sub").get_to(s.sub_indices);
    if (j.contains("ref")) {
        j.at("ref").get_to(s.ref_indices);
    }
    if (j.contains("reverse")) {
        j.at("reverse").get_to(s.reverse);
    }
    if (j.contains("inverse")) {
        j.at("inverse").get_to(s.inverse);
    }
    if (j.contains("scale")) {
        j.at("scale").get_to(s.scale);
    }
}

void to_json(json &j, const MTContrast &s) {
    j = json{
        {"name", s.name},
        {"add", s.add_indices},
        {"sub", s.sub_indices},
        {"ref", s.ref_indices},
        {"reverse", s.reverse},
    };
}

/*
 * Main
 */
int mtr_main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser(
        "Calculates MTR/ihMTR/MTasym etc. By default calculate MTR assuming input file has two"
        "volumes, first MTw second Ref\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(
        parser, "INPUT", "Input file with different MT contrasts");

    args::HelpFlag       help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=hardware limit or $QUIT_THREADS)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(
        parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    args::ValueFlag<std::string> mask_path(
        parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> json_file(
        parser, "JSON", "Read custom contrasts from JSON file", {"json"});
    args::ValueFlag<std::string> reference(parser, "REF", "External reference image", {"ref", 'r'});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    auto const input_img = QI::ReadImage<QI::VectorVolumeF>(QI::CheckPos(input_path), verbose);
    QI::VolumeF::Pointer const ref_img =
        reference ? QI::ReadImage(reference.Get(), verbose) : nullptr;

    std::vector<MTContrast> contrasts;
    if (json_file) {
        QI::Log(verbose, "Reading contrasts");
        json doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
        doc["contrasts"].get_to(contrasts);
    } else {
        contrasts.push_back({"MTR", {0}, {}, {1}, true});
    }

    QI::VolumeF::Pointer mask_img = nullptr;
    if (mask_path) {
        mask_img = QI::ReadImage(mask_path.Get(), verbose);
    }

    QI::Info(verbose, "Allocating output memory");
    std::vector<QI::VolumeF::Pointer> out_imgs;
    for (size_t ii = 0; ii < contrasts.size(); ii++) {
        out_imgs.push_back(QI::NewImageLike<QI::VolumeF>(input_img));
    }

    QI::Info(verbose, "Processing");
    auto mt = itk::MultiThreaderBase::New();
    mt->SetNumberOfWorkUnits(threads.Get());
    mt->ParallelizeImageRegion<3>(
        input_img->GetBufferedRegion(),
        [&](const QI::VolumeF::RegionType &region) {
            itk::ImageRegionConstIterator<QI::VectorVolumeF> in_it(input_img, region);
            itk::ImageRegionConstIterator<QI::VolumeF>       ref_it;
            if (ref_img) {
                ref_it = itk::ImageRegionConstIterator<QI::VolumeF>(ref_img, region);
            }
            std::vector<itk::ImageRegionIterator<QI::VolumeF>> out_its;
            for (auto const &img : out_imgs) {
                out_its.push_back(itk::ImageRegionIterator<QI::VolumeF>(img, region));
            }
            itk::ImageRegionConstIterator<QI::VolumeF> mask_it;
            if (mask_img) {
                mask_it = itk::ImageRegionConstIterator<QI::VolumeF>(mask_img, region);
            }

            for (in_it.GoToBegin(); !in_it.IsAtEnd(); ++in_it) {
                if (!mask_img || mask_it.Get()) {
                    const Eigen::Map<const Eigen::ArrayXf> in(
                        in_it.Get().GetDataPointer(), input_img->GetNumberOfComponentsPerPixel());
                    for (size_t ii = 0; ii < out_its.size(); ii++) {
                        auto &it  = out_its[ii];
                        auto &con = contrasts[ii];
                        float pos = 0.f, neg = 0.f, ref = 0.f;
                        for (auto const &ind : con.add_indices) {
                            auto const &pos_val = in[ind];
                            pos += con.inverse ? 1. / (pos_val * con.add_indices.size()) :
                                                 pos_val / con.add_indices.size();
                        }
                        for (auto const &ind : con.sub_indices) {
                            auto const &neg_val = in[ind];
                            neg += con.inverse ? 1. / (neg_val * con.sub_indices.size()) :
                                                 neg_val / con.sub_indices.size();
                        }

                        for (auto const &ind : con.ref_indices) {
                            if (ind == -1 && ref_img) {
                                ref += ref_it.Get() / con.ref_indices.size();
                            } else {
                                ref += in[ind] / con.ref_indices.size();
                            }
                        }
                        float const diff    = con.scale * (pos - neg);
                        float const value   = con.reverse ? ref - diff : diff;
                        float const percent = 100.f * (con.inverse ? (ref * value) : (value / ref));
                        it.Set(percent);
                    }
                } else {
                    for (auto &it : out_its) {
                        it.Set(0.0);
                    }
                }
                for (auto &it : out_its) {
                    ++it;
                }
                if (mask_img) {
                    ++mask_it;
                }
                if (ref_img) {
                    ++ref_it;
                }
            }
        },
        nullptr);

    QI::Info(verbose, "Finished");
    for (size_t ii = 0; ii < contrasts.size(); ii++) {
        QI::WriteImage(out_imgs[ii], outarg.Get() + contrasts[ii].name + QI::OutExt(), verbose);
    }
    return EXIT_SUCCESS;
}
