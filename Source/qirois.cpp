/*
 *  qirois.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <fstream>
#include <string>

#include "QI/Types.h"
#include "QI/Util.h"
#include "QI/Option.h"

#include "itkLabelGeometryImageFilter.h"
typedef itk::LabelGeometryImageFilter<QI::VolumeI, QI::VolumeF> TLblGeoFilter;

const std::string usage{
"Usage is: qirois [options] labelfile\n\
\n\
Calculates volumes of all ROI labels in a file. Uses itkLabelGeometryImageFilter\n"
};

int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::OptionList opts(usage);
    QI::Option<std::string> label_list("",'l',"labels","Specify labels to search for", opts);
    QI::Switch print_labels('p',"print_labels","Print out label/ROI numbers first", opts);
    QI::Help help(opts);
    std::vector<std::string> nonopts = opts.parse(argc, argv);
    if (nonopts.size() != 1) {
        std::cerr << opts << std::endl;
        std::cerr << "No input filename specified." << std::endl;
        return EXIT_FAILURE;
    }
    QI::VolumeI::Pointer label_image = QI::ReadImage<QI::VolumeI>(nonopts[0]);
    double voxvol = label_image->GetSpacing()[0]; for (int i = 1; i < 3; i++) voxvol *= label_image->GetSpacing()[i];
    TLblGeoFilter::Pointer label_filter = TLblGeoFilter::New();
    label_filter->SetInput(label_image);
    label_filter->CalculatePixelIndicesOff();
    label_filter->CalculateOrientedBoundingBoxOff();
    label_filter->CalculateOrientedLabelRegionsOff();
    label_filter->Update();

    typename TLblGeoFilter::LabelsType labels;
    if (label_list.set()) {
        std::ifstream file(*label_list);
        std::string line;
        while (std::getline(file, line)) {
            std::stringstream linestream(line);
            int label;
            linestream >> label;
            labels.push_back(label);
        }
    } else {
        labels = label_filter->GetLabels();
        std::sort(labels.begin(), labels.end());
    }
    for(auto it = labels.begin(); it != labels.end(); it++) {
        if (*print_labels) std::cout << *it << ",";
        std::cout << (voxvol * label_filter->GetVolume(*it)) << std::endl;
    }
    return EXIT_SUCCESS;
}

