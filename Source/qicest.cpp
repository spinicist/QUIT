/*
 *  qicest.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "QI/Util.h"
#include "QI/Option.h"
#include "QI/CEST.h"
#include "Filters/ApplyAlgorithmFilter.h"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::OptionList opts("Usage is: qicest [options] input_file\n\nA utility for calculating CEST metrics.\nInput must be a single image.");
    QI::Option<int> num_threads(4,'T',"threads","Use N threads (default=4, 0=hardware limit)", opts);
    QI::RegionOption<QI::ApplyF::TRegion> subregion('s',"subregion","Process subregion starting at voxel I,J,K with size SI,SJ,SK", opts);
    QI::ImageOption<QI::VolumeF> mask('m', "mask", "Mask input with specified file", opts);
    QI::Option<std::string> outPrefix("", 'o', "out","Prefix output filenames", opts);
    QI::Switch suppress('n',"no-prompt","Suppress input prompts", opts);
    QI::Switch verbose('v',"verbose","Print more information", opts);
    QI::Help help(opts);
    std::deque<std::string> nonopts = opts.parse(argc, argv);
    if (nonopts.size() != 1) {
        std::cerr << opts << std::endl;
        std::cerr << "No input filename specified." << std::endl;
        return EXIT_FAILURE;
    }
    if (*verbose) cout << "Opening file: " << nonopts[0] << endl;
    auto data = QI::ReadVectorImage<float>(nonopts[0]);

    if (!*suppress) cout << "Enter input frequencies: " << endl;
    Eigen::ArrayXf ifrqs; QI::ReadArray(cin, ifrqs);
    if (!*suppress) cout << "Enter output spectrum frequencies: " << endl;
    Eigen::ArrayXf ofrqs; QI::ReadArray(cin, ofrqs);
    if (!*suppress) cout << "Enter output asymmetry frequencies: " << endl;
    Eigen::ArrayXf afrqs; QI::ReadArray(cin, afrqs); // Asymmetry output
    shared_ptr<QI::CESTAlgo> algo = make_shared<QI::CESTAlgo>(ifrqs, ofrqs, afrqs);

    auto apply = QI::ApplyVectorF::New();
    apply->SetAlgorithm(algo);
    apply->SetPoolsize(*num_threads);
    apply->SetInput(0, data);
    apply->SetMask(*mask);
    if (subregion.set()) {
        if (*verbose) cout << "Setting subregion: " << *subregion << endl;
        apply->SetSubregion(*subregion);
    }
    if (*verbose) {
        cout << "Processing" << endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (*verbose) {
        cout << "Elapsed time was " << apply->GetTotalTime() << "s" << endl;
        cout << "Writing results files." << endl;
    }
    *outPrefix = *outPrefix + "CEST_";
    for (int i = 0; i < algo->numOutputs(); i++) {
        QI::WriteVectorImage(apply->GetOutput(i), *outPrefix + algo->names().at(i) + QI::OutExt());
    }
    
    if (*verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
