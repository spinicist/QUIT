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
#include "itkMultiThreaderBase.h"

/*
 * Main
 */
int noise_est_main(args::Subparser &parser) {
    args::Positional<std::string> input_path(parser, "FILE", "Input 4D file");

    args::ValueFlag<std::string> region(
        parser, "REGION", "Measure noise in specified region", {'r', "region"});
    args::ValueFlag<std::string> mask(
        parser, "MASK", "Measure noise in specified mask", {'m', "mask"});
    args::Flag mean_sqr(parser,
                        "MEAN^2",
                        "Return the mean of the squared values (for Rician correction)",
                        {"meansqr"});
    parser.Parse();
    if (!region && !mask) {
        QI::Fail("Either REGION or MASK must be specified");
    }

    auto const input  = QI::ReadImage<QI::VectorVolumeF>(QI::CheckPos(input_path), verbose);
    int const  insize = input->GetNumberOfComponentsPerPixel();
    QI::VolumeF::Pointer mask_img = mask ? QI::ReadImage(mask.Get(), verbose) : nullptr;

    double noise_mean = 0., noise_sqr_mean = 0., noise_sigma = 0.;

    // As per https://linkinghub.elsevier.com/retrieve/pii/0730725X93902253
    QI::Log(verbose, "Calculating noise statistics");
    auto mt           = itk::MultiThreaderBase::New();
    auto noise_region = region ? QI::RegionFromString<QI::VectorVolumeF::RegionType>(region.Get()) :
                                 input->GetLargestPossibleRegion();
    itk::ImageRegionConstIterator<QI::VectorVolumeF> noise_it(input, noise_region);
    itk::ImageRegionConstIterator<QI::VolumeF>       mask_it;
    if (mask) {
        mask_it = itk::ImageRegionConstIterator<QI::VolumeF>(mask_img, noise_region);
    }
    long samples = 0;
    for (noise_it.GoToBegin(); !noise_it.IsAtEnd(); ++noise_it) {
        if (!mask || mask_it.Get()) {
            Eigen::Map<Eigen::VectorXf const> const noise_vec(noise_it.Get().GetDataPointer(),
                                                              insize);
            for (Eigen::Index ii = 0; ii < insize; ii++) {
                auto const val = noise_vec[ii];
                noise_mean += val;
                noise_sqr_mean += val * val;
                samples++;
            }
        }
        if (mask) {
            ++mask_it;
        }
    }
    noise_mean     = noise_mean / samples;
    noise_sqr_mean = noise_sqr_mean / samples;
    noise_sigma    = sqrt(noise_sqr_mean - noise_mean * noise_mean);
    QI::Log(verbose,
            "Noise samples {} mean {:.3g} sd {:.3g} ratio {:.3g} squared mean {:.3g}",
            samples,
            noise_mean,
            noise_sigma,
            noise_mean / noise_sigma,
            noise_sqr_mean);

    if (mean_sqr) {
        fmt::print("{}\n", noise_sqr_mean);
    } else {
        fmt::print("{}\n", noise_sigma);
    }

    return EXIT_SUCCESS;
}
