/*
 *  qiafi.cpp
 *
 *  Created by Tobias Wood on 2015/06/04.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <string>

#include "Args.h"
#include "ImageIO.h"
#include "ImageTypes.h"
#include "Util.h"

#include "itkBinaryFunctorImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkExtractImageFilter.h"

template <class TPixel> class AFI {
  public:
    AFI() {}
    ~AFI() {}
    bool operator!=(const AFI &) const { return false; }
    bool operator==(const AFI &other) const { return !(*this != other); }

    inline TPixel operator()(const TPixel &r, const TPixel &n) const {
        TPixel temp = (r * n - 1.) / (n - r);
        if (temp > 1.)
            temp = 1.;
        if (temp < -1.)
            temp = -1.;
        TPixel alpha = acos(temp) * 180. / M_PI;
        return alpha;
    }
};

int afi_main(int argc, char **argv) {
    args::ArgumentParser parser(
        "Calculates B1 maps from AFI data. Input file should have two volumes\n"
        "http://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "INPUT", "Input file");
    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=hardware limit or $QUIT_THREADS)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string> out_prefix(
        parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<double> nom_flip(
        parser, "NOMINAL FLIP", "Specify nominal flip-angle, default 55", {'f', "flip"}, 55.0);
    args::ValueFlag<double> tr_ratio(
        parser, "TR RATIO", "Specify TR2:TR1 ratio, default 5", {'r', "ratio"}, 5.0);
    args::Flag save_angle(
        parser, "SAVE ANGLE", "Write out the actual flip-angle as well as B1", {'s', "save"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    auto inFile = QI::ReadImage<QI::SeriesF>(QI::CheckPos(input_path), verbose);
    QI::Log(verbose, "Nominal flip-angle = {} degrees", nom_flip.Get());
    QI::Log(verbose, "TR2:TR1 ratio = {}", tr_ratio.Get());
    auto volume1                   = itk::ExtractImageFilter<QI::SeriesF, QI::VolumeF>::New();
    auto volume2                   = itk::ExtractImageFilter<QI::SeriesF, QI::VolumeF>::New();
    auto region                    = inFile->GetLargestPossibleRegion();
    region.GetModifiableSize()[3]  = 0;
    region.GetModifiableIndex()[3] = 0;
    volume1->SetExtractionRegion(region);
    volume1->SetInput(inFile);
    volume1->SetDirectionCollapseToSubmatrix();
    region.GetModifiableIndex()[3] = 1;
    volume2->SetExtractionRegion(region);
    volume2->SetInput(inFile);
    volume2->SetDirectionCollapseToSubmatrix();

    auto imageRatio = itk::DivideImageFilter<QI::VolumeF, QI::VolumeF, QI::VolumeF>::New();
    imageRatio->SetInput(0, volume2->GetOutput());
    imageRatio->SetInput(1, volume1->GetOutput());
    auto afi =
        itk::BinaryFunctorImageFilter<QI::VolumeF, QI::VolumeF, QI::VolumeF, AFI<float>>::New();
    afi->SetInput1(imageRatio->GetOutput());
    afi->SetConstant2(tr_ratio.Get());
    auto B1 = itk::DivideImageFilter<QI::VolumeF, QI::VolumeF, QI::VolumeF>::New();
    B1->SetInput1(afi->GetOutput());
    B1->SetConstant2(nom_flip.Get());
    B1->Update();
    QI::WriteImage(B1->GetOutput(), out_prefix.Get() + "AFI_B1" + QI::OutExt(), verbose);
    if (save_angle)
        QI::WriteImage(afi->GetOutput(), out_prefix.Get() + "AFI_angle" + QI::OutExt(), verbose);
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
