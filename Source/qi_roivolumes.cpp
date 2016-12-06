/*
 *  qi_roivolumes.cpp
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
"Usage is: qi_roivolumes [options] label_image1 [ ... label_imageN] \n\
\n\
Calculates volumes of ROI labels in files. Uses itkLabelGeometryImageFilter\n"
};

double VoxelVolume(QI::VolumeI::Pointer img) {
    double vox_volume = img->GetSpacing()[0];
    for (int v = 1; v < 3; v++)
        vox_volume *= img->GetSpacing()[v];
    return vox_volume;
}

int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::OptionList opts(usage);
    QI::Option<std::string> label_list("",'l',"labels","Specify labels to use in text file", opts);
    QI::Switch ignore_zero('z',"ignore_zero","Ignore 0 label (background)", opts);
    QI::Switch print_labels('p',"print_labels","Print out label/ROI numbers first", opts);
    QI::Switch rows('r',"rows","Print comma-separated rows instead of columns", opts);
    QI::Help help(opts);
    std::deque<std::string> nonopts = opts.parse(argc, argv);
    if (nonopts.size() < 1) {
        std::cerr << opts << std::endl;
        std::cerr << "No input filename specified." << std::endl;
        return EXIT_FAILURE;
    }
    
    TLblGeoFilter::Pointer filter = TLblGeoFilter::New();
    filter->CalculatePixelIndicesOff();
    filter->CalculateOrientedBoundingBoxOff();
    filter->CalculateOrientedLabelRegionsOff();
    QI::VolumeI::Pointer img = QI::ReadImage<QI::VolumeI>(nonopts.at(0));
    filter->SetInput(img);
    filter->Update();
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
        labels = filter->GetLabels();
        std::sort(labels.begin(), labels.end());
    }
    if (*ignore_zero) {
        std::remove(labels.begin(), labels.end(), 0);
    }

    double vox_volume = VoxelVolume(img);
    std::vector<std::vector<double>> volumes(nonopts.size(), std::vector<double>(labels.size()));

    for (int i = 0; i < labels.size(); i++) {
        volumes.at(0).at(i) = vox_volume * filter->GetVolume(labels.at(i));
    }

    for (int f = 1; f < nonopts.size(); f++) {
        img = QI::ReadImage<QI::VolumeI>(nonopts.at(f));
        filter->SetInput(img);
        filter->Update();
        vox_volume = VoxelVolume(img);
        for (int i = 0; i < labels.size(); i++) {
            volumes.at(f).at(i) = vox_volume * filter->GetVolume(labels.at(i));
        }
    }

    if (*rows) {
        if (*print_labels) {
            auto it = labels.begin();
            std::cout << *it;
            for (++it; it != labels.end(); ++it) {
                std::cout << "," << *it;
            }
            std::cout << std::endl;
        }
        for (auto f = volumes.begin(); f != volumes.end(); ++f){
            auto v = f->begin();
            std::cout << *v;
            for (++v; v != f->end(); ++v) {
                std::cout << "," << *v;
            }
            std::cout << std::endl;
        }
    } else {
        for (int l = 0; l < labels.size(); ++l) {
            if (*print_labels)
                std::cout << labels.at(l) << ",";
            
            auto f = volumes.begin();
            std::cout << f->at(l);
            for (++f; f != volumes.end(); ++f) {
                std::cout << "," << f->at(l);
            }
            std::cout << std::endl;
        }
    }
    return EXIT_SUCCESS;
}

