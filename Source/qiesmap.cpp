/*
 *  qiesmap.cpp
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
#include "QI/Algorithms/Banding.h"
#include "QI/Algorithms/Ellipse.h"
#include "QI/Algorithms/EllipseFit.h"
#include "Filters/ApplyAlgorithmFilter.h"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::OptionList opts("Usage is: qiesmap [options] input_file\n\nA utility for calculating T1,T2,PD and f0 maps from SSFP data.\nInput must be a single complex image with at least 6 phase increments.");
    QI::Option<int> num_threads(4,'T',"threads","Use N threads (default=4, 0=hardware limit)", opts);
    QI::RegionOption<QI::ApplyF::TRegion> subregion('s',"subregion","Process subregion starting at voxel I,J,K with size SI,SJ,SK", opts);
    QI::ImageOption<QI::VolumeF> mask('m', "mask", "Mask input with specified file", opts);
    QI::ImageOption<QI::VolumeF> B1('b', "B1", "B1 Map file (ratio)", opts);
    QI::Option<int> ph_incs(6,'\0',"ph_incs","Number of phase increments (default is 6).", opts);
    QI::Switch ph_order('\0',"ph_order","Data order is phase, then flip-angle (default opposite).", opts);
    QI::Switch alt_order('\0',"alt_order","Opposing phase-incs alternate (default is 2 blocks)", opts);
    QI::EnumOption algorithm("hcf",'h','a',"algo","Choose algorithm (h/c/f)", opts);
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
    auto data = QI::ReadVectorImage<complex<float>>(nonopts[0]);
    shared_ptr<QI::ESAlgo> algo;
    auto seq = make_shared<QI::SSFPEcho>(cin, !*suppress);
    if (*verbose) cout << *seq;
    switch (*algorithm) {
    case 'h': algo = make_shared<QI::HyperEllipse>(seq, *ph_order); break;
    case 'c': algo = make_shared<QI::ConstrainedEllipse>(seq, *ph_order, *alt_order); break;
    case 'f': algo = make_shared<QI::FitEllipse>(seq, *ph_order, *alt_order); break;
    }
    auto apply = QI::ApplyVectorXFVectorF::New();
    apply->SetAlgorithm(algo);
    apply->SetPoolsize(*num_threads);
    apply->SetSplitsPerThread(*num_threads);
    apply->SetInput(0, data);
    apply->SetMask(*mask);
    apply->SetConst(0, *B1);
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
    *outPrefix = *outPrefix + "ES_";
    for (int i = 0; i < algo->numOutputs(); i++) {
        QI::WriteVectorImage(apply->GetOutput(i), *outPrefix + algo->names().at(i) + QI::OutExt());
    }
    
    if (*verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
