/*
 *  qi_pca.cpp
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>

#include <Eigen/Eigenvalues>

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkMultiThreaderBase.h"

// #define QI_DEBUG_BUILD 1

#include "Args.h"
#include "ImageIO.h"
#include "JSON.h"
#include "Macro.h"
#include "Spline.h"
#include "Util.h"

using namespace std::literals;

int pca_main(args::Subparser &parser) {
    args::Positional<std::string> input_path(parser, "INPUT", "Input 4D file");
    args::ValueFlag<int>          threads(parser,
                                 "THREADS",
                                 "Use N threads (default=hardware limit or $QUIT_THREADS)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string>  outarg(
        parser, "OUTPUT", "Change ouput filename (default is input_pca)", {'o', "out"});
    args::ValueFlag<std::string> project(
        parser, "PROJECT", "Save proj_img onto PCs to specified 4D file", {'p', "project"});
    args::ValueFlag<std::string> save_pcs(
        parser, "PC JSON", "Save PCs into specified JSON file", {'s', "save_pcs"});
    args::ValueFlag<int> n_retain(
        parser, "RETAIN", "Number of PCs to retain, default 3", {'r', "retain"}, 3);
    args::ValueFlag<std::string> mask(
        parser, "MASK", "Only process voxels within the mask (recommended)", {'m', "mask"});
    parser.Parse();

    auto const input  = QI::ReadImage<QI::VectorVolumeF>(QI::CheckPos(input_path), verbose);
    auto const region = input->GetLargestPossibleRegion();

    QI::VolumeF::Pointer const mask_img = mask ? QI::ReadImage(mask.Get(), verbose) : nullptr;
    auto const                 Nvox     = [&]() {
        if (mask) {
            QI::Log(verbose, "Counting voxels in mask...");
            auto         mask_it = itk::ImageRegionConstIterator<QI::VolumeF>(mask_img, region);
            Eigen::Index Nv = 0;
            while (!mask_it.IsAtEnd()) {
                if (mask_it.Get())
                    ++Nv;
                ++mask_it;
            }
            return Nv;
        } else {
            QI::Log(verbose, "No mask, will use all voxels in image");
            return static_cast<Eigen::Index>(region.GetNumberOfPixels());
        }
    }();
    QI::Log(verbose, "Total voxels = {}", Nvox);

    Eigen::Index const Nq   = input->GetNumberOfComponentsPerPixel();
    Eigen::Index const Nret = (n_retain.Get() > Nq) ? Nq : n_retain.Get();

    Eigen::MatrixXd X(Nvox, Nq);
    {
        // Use scope to keep these iterators separate to the multi-threaded ones later
        itk::ImageRegionConstIterator<QI::VectorVolumeF> in_it(input, region);
        itk::ImageRegionConstIterator<QI::VolumeF>       mask_it;
        if (mask)
            mask_it = itk::ImageRegionConstIterator<QI::VolumeF>(mask_img, region);
        Eigen::Index row_ind = 0;
        for (in_it.GoToBegin(); !in_it.IsAtEnd(); ++in_it) {
            if (!mask || mask_it.Get()) {
                const Eigen::Map<const Eigen::VectorXf> vox(in_it.Get().GetDataPointer(), Nq);
                X.row(row_ind) = vox.cast<double>();
                ++row_ind;
            }
            if (mask) {
                ++mask_it;
            }
        }
    }

    QI::Info(verbose, "Calculating Principal Components");
    Eigen::VectorXd xmean = X.colwise().mean();
    X.rowwise() -= xmean.transpose();
    auto const cov = X.adjoint() * X / (Nq - 1);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
    auto const norm_eig_vals = eig.eigenvalues() / eig.eigenvalues().sum();

    Eigen::VectorXd retained_vals = norm_eig_vals.tail(Nret).reverse();
    Eigen::MatrixXd retained_vecs =
        eig.eigenvectors().rightCols(Nret).rowwise().reverse().colwise().normalized();
    QI::Log(verbose,
            "Retaining {} eigenvalues with % variance: {}",
            Nret,
            (retained_vals * 100).transpose());

    if (save_pcs) {
        json doc;
        doc["eigenvalues"] = retained_vals;
        for (Eigen::Index v = 0; v < Nret; v++) {
            // std::string name = "eigenvector_"s + std::to_string(v);
            doc["eigenvector"] = retained_vecs.col(v);
        }
        QI::Log(verbose, "Saving PCs to JSON file: {}", save_pcs.Get());
        QI::WriteJSON(save_pcs.Get(), doc);
    }

    auto proj_img = QI::NewImageLike<QI::VectorVolumeF>(input, Nret);
    auto out_img  = QI::NewImageLike<QI::VectorVolumeF>(input, Nq);

    QI::Info(verbose, "Calculating projection...");
    auto mt = itk::MultiThreaderBase::New();
    mt->SetNumberOfWorkUnits(threads.Get());
    mt->ParallelizeImageRegion<3>(
        input->GetBufferedRegion(),
        [&](const QI::VectorVolumeF::RegionType &thread_region) {
            itk::ImageRegionConstIterator<QI::VectorVolumeF> in_it(input, thread_region);
            itk::ImageRegionIterator<QI::VectorVolumeF>      proj_it(proj_img, thread_region);
            itk::ImageRegionIterator<QI::VectorVolumeF>      out_it(out_img, thread_region);
            itk::ImageRegionConstIterator<QI::VolumeF>       mask_it;
            if (mask_img)
                mask_it = itk::ImageRegionConstIterator<QI::VolumeF>(mask_img, thread_region);
            for (in_it.GoToBegin(); !in_it.IsAtEnd(); ++in_it, ++proj_it, ++out_it) {
                itk::VariableLengthVector<float> itk_proj(Nret);
                itk::VariableLengthVector<float> itk_out(Nq);
                if (!mask_img || mask_it.Get()) {
                    const Eigen::Map<const Eigen::VectorXf> vox(in_it.Get().GetDataPointer(), Nq);

                    auto const &    voxd     = vox.cast<double>() - xmean;
                    Eigen::VectorXd vox_proj = voxd.transpose() * retained_vecs;
                    Eigen::VectorXd vox_out  = xmean + retained_vecs * vox_proj;

                    for (Eigen::Index i = 0; i < Nret; i++) {
                        itk_proj[i] = vox_proj[i];
                    }
                    for (Eigen::Index i = 0; i < Nq; i++) {
                        itk_out[i] = vox_out[i];
                    }
                } else {
                    itk_proj.Fill(0.0);
                    itk_out.Fill(0.0);
                }
                proj_it.Set(itk_proj);
                out_it.Set(itk_out);

                if (mask) {
                    ++mask_it;
                }
            }
        },
        nullptr);
    QI::Info(verbose, "Finished");
    if (project) {
        QI::WriteImage(proj_img, project.Get(), verbose);
    }
    std::string outname = outarg ?
                              outarg.Get() :
                              QI::StripExt(QI::Basename(input_path.Get())) + "_pca" + QI::OutExt();
    QI::WriteImage(out_img, outname, verbose);
    return EXIT_SUCCESS;
}