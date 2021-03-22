/*
 *  qi_denoise.cpp
 *
 *  Copyright (c) 2020 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "tgv-denoise.hpp"

#include "Args.h"
#include "ImageIO.h"
#include "Util.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImportImageFilter.h"
#include "itkMultiThreaderBase.h"

/*
 * Main
 */
int denoise_main(args::Subparser &parser) {
    args::Positional<std::string> iname(parser, "INPUT", "Input filename");
    args::Positional<std::string> oname(parser, "OUTPUT", "Output filename");

    args::ValueFlag<long> its(
        parser, "MAX ITS", "Maximum number of iterations (16)", {'i', "max_its"}, 16);
    args::ValueFlag<float> thr(
        parser, "TRESHOLD", "Threshold for termination (1e-10)", {"thresh"}, 1.e-10);
    args::ValueFlag<float> alpha(
        parser, "ALPHA", "Regularisation weighting (1e-5)", {"alpha", 'a'}, 1.e-5f);
    args::ValueFlag<float> alpha_reduction(
        parser, "ALPHA REDUCE", "Reduce by factor (suggest 0.1)", {"reduce", 'r'}, 1.f);
    args::ValueFlag<float> step_size(
        parser, "STEP SIZE", "Inverse of step size (default 8)", {"step"}, 8.f);
    args::Flag complex(parser, "COMPLEX", "Input is complex-valued", {"complex", 'x'});
    parser.Parse();
    if (!iname && !oname) {
        QI::Fail("Both input and output filename must be set");
    }

    Eigen::ThreadPool       pool(std::thread::hardware_concurrency());
    Eigen::ThreadPoolDevice device(&pool, pool.NumThreads());

    auto pipeline = [&]<typename T>() {
        using TT        = Eigen::Tensor<T, 4>;
        auto const iimg = QI::ReadImage<itk::Image<T, 4>>(iname.Get(), verbose);
        auto const sz   = iimg->GetLargestPossibleRegion().GetSize();

        typename TT::Dimensions dims;
        std::copy_n(sz.begin(), 4, dims.begin());
        Eigen::TensorMap<TT> input(iimg->GetBufferPointer(), dims);
        TT                   output(dims);
        for (long ii = 0; ii < dims[3]; ii++) {
            QI::Log(verbose, "Processing volume {}", ii);
            Eigen::Tensor<T, 3> vol     = input.template chip<3>(ii);
            output.template chip<3>(ii) = tgvdenoise(vol,
                                                     its.Get(),
                                                     thr.Get(),
                                                     alpha.Get(),
                                                     alpha_reduction.Get(),
                                                     step_size.Get(),
                                                     verbose,
                                                     device);
        }
        using Importer                    = itk::ImportImageFilter<T, 4>;
        typename Importer::Pointer import = Importer::New();
        import->SetRegion(iimg->GetLargestPossibleRegion());
        import->SetSpacing(iimg->GetSpacing());
        import->SetOrigin(iimg->GetOrigin());
        import->SetDirection(iimg->GetDirection());
        import->SetImportPointer(output.data(), output.size(), false);
        import->Update();
        QI::WriteImage(import->GetOutput(), oname.Get(), verbose);
    };

    if (complex) {
        pipeline.operator()<std::complex<float>>();
    } else {
        pipeline.operator()<float>();
    }

    return EXIT_SUCCESS;
}
