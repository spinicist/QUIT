/*
 *  qidream.cpp
 *
 *  Created by Tobias Wood on 2015/06/04.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <string>
#include <getopt.h>

#include "QI/Types.h"
#include "QI/Util.h"
#include "QI/IO.h"

#include "itkExtractImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"

template<class TPixel> class DREAM {
public:
    DREAM() {}
    ~DREAM() {}
    bool operator!=(const DREAM &) const { return false; }
    bool operator==(const DREAM &other) const { return !(*this != other); }

    inline TPixel operator()(const TPixel &p1,
                             const TPixel &p2) const
    {
        TPixel alpha = atan(sqrt(2.*p1/p2)) * 180. / M_PI;
        return alpha;
    }
};

int main(int argc, char **argv) {
    QI::ArgParser args{argc, argv,
        "Usage is: qidream input [options]\n"
        "Calculates a B1 (flip-angle) map from DREAM data. Input must have 2 volumes.",
        {{"help",      'h', "Display the help message and quit", false},
         {"verbose",   'v', "Print more information", false},
         {"no-prompt", 'n', "Suppress input prompts", false},
         {"threads",   'T', "Use N threads (default=4, 0=hardware limit)", true},
         {"out",       'o', "Add a prefix to output filenames", true},
         {"order",     'O', "Volume order - f/s/v - fid/ste/vst first", true},
         {"alpha",     'a', "Reference flip-angle (default 55)", true},
         {"subregion", 's', "Process subregion starting at voxel I,J,K with size SI,SJ,SK", false}}
    };

    bool verbose = args.option_present("verbose");
    std::string outPrefix = args.string_value("out","") + "DREAM_";
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(args.option_value("threads",4));
    if (args.nonoptions().size() != 1) {
        std::cerr << "Must specify one input file with two volumes." << std::endl;
        return EXIT_FAILURE;
    }
    if (verbose) std::cout << "Opening input file " << args.nonoptions()[0] << std::endl;
    auto inFile = QI::ReadImage<QI::SeriesF>(args.nonoptions()[0]);

    auto volume1 = itk::ExtractImageFilter<QI::SeriesF, QI::VolumeF>::New();
    auto volume2 = itk::ExtractImageFilter<QI::SeriesF, QI::VolumeF>::New();
    auto region = inFile->GetLargestPossibleRegion();
    region.GetModifiableSize()[3] = 0;
    switch (args.option_value("order", 'f')) {
    case 'f': region.GetModifiableIndex()[3] = 0; break;
    case 's':
    case 'v': region.GetModifiableIndex()[3] = 1; break;
    }
    volume1->SetExtractionRegion(region);
    volume1->SetInput(inFile);
    volume1->SetDirectionCollapseToSubmatrix();
    // Swap to other volume
    region.GetModifiableIndex()[3] = (region.GetIndex()[3] + 1) % 2;
    volume2->SetExtractionRegion(region);
    volume2->SetInput(inFile);
    volume2->SetDirectionCollapseToSubmatrix();

    auto dream = itk::BinaryFunctorImageFilter<QI::VolumeF,
                                               QI::VolumeF,
                                               QI::VolumeF,
                                               DREAM<float>>::New();
    dream->SetInput1(volume1->GetOutput());
    dream->SetInput2(volume2->GetOutput());

    auto B1 = itk::DivideImageFilter<QI::VolumeF, QI::VolumeF, QI::VolumeF>::New();
    B1->SetInput1(dream->GetOutput());
    B1->SetConstant2(args.option_value("alpha", 55));
    B1->Update();
    QI::WriteImage(dream->GetOutput(), outPrefix + "angle" + QI::OutExt());
    QI::WriteImage(B1->GetOutput(),  outPrefix + "B1" + QI::OutExt());
    if (verbose) std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}

