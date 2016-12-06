/*
 *  qi_roivalues.cpp
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
"Usage is: qi_roivalues [options] label_image val_image1 [ ... val_imageN] \n\
\n\
Extracts average values of ROI labels in image files. Uses itkLabelGeometryImageFilter\n\
\n\
ROIs with a volume of 0 will have NaN average values.\n\
Missing files will issue a warning and insert NaN for all ROIs\n"
};

double VoxelVolume(QI::VolumeI::Pointer lbl_img) {
    double vox_volume = lbl_img->GetSpacing()[0];
    for (int v = 1; v < 3; v++)
        vox_volume *= lbl_img->GetSpacing()[v];
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
    std::deque<std::string> image_list = opts.parse(argc, argv);
    if (image_list.size() < 2) {
        std::cerr << opts << std::endl;
        std::cerr << "Must specify label image and input image." << std::endl;
        return EXIT_FAILURE;
    }
    
    TLblGeoFilter::Pointer filter = TLblGeoFilter::New();
    filter->CalculatePixelIndicesOff();
    filter->CalculateOrientedBoundingBoxOff();
    filter->CalculateOrientedLabelRegionsOff();
    QI::VolumeI::Pointer lbl_img = QI::ReadImage<QI::VolumeI>(image_list.at(0));
    image_list.pop_front();
    filter->SetInput(lbl_img);

    // Update once to get label list
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

    std::vector<std::vector<double>> mean_values(image_list.size(), std::vector<double>(labels.size()));

    // Now loop through value images and get average 'intensity' for each ROI
    for (int i = 0; i < image_list.size(); i++) {
        QI::VolumeF::Pointer val_img;
        try {
            val_img = QI::ReadImage(image_list.at(i));
            filter->SetIntensityInput(val_img);
            filter->Update();
            for (int l = 0; l < labels.size(); l++) {
                const double value = filter->GetIntegratedIntensity(labels.at(l));
                const double volume = filter->GetVolume(labels.at(l));
                mean_values.at(i).at(l) =  value / volume;
            }
        } catch (itk::ExceptionObject &err) {
            std::cerr << "Failed to read file: " << image_list.at(i) << std::endl;
            std::fill(mean_values.at(i).begin(), mean_values.at(i).end(), std::nan(""));
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
        for (auto f = mean_values.begin(); f != mean_values.end(); ++f){
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
            
            auto f = mean_values.begin();
            std::cout << f->at(l);
            for (++f; f != mean_values.end(); ++f) {
                std::cout << "," << f->at(l);
            }
            std::cout << std::endl;
        }
    }
    return EXIT_SUCCESS;
}

