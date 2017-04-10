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
#include "QI/IO.h"
#include "QI/Args.h"

#include "itkLabelGeometryImageFilter.h"
typedef itk::LabelGeometryImageFilter<QI::VolumeI, QI::VolumeF> TLblGeoFilter;

double VoxelVolume(QI::VolumeI::Pointer img) {
    double vox_volume = img->GetSpacing()[0];
    for (int v = 1; v < 3; v++)
        vox_volume *= img->GetSpacing()[v];
    return vox_volume;
}

int main(int argc, char **argv) {
    Eigen::initParallel();

    args::ArgumentParser parser("Calculates volumes of ROI labels. Multiple inputs can be specified.\n"
                                "Output goes to stdout\n"
                                "http://github.com/spinicist/QUIT");
    args::PositionalList<std::string> in_paths(parser, "INPUT", "Input file paths.");
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information (will mess up output, for test runs only)", {'v', "verbose"});
    args::ValueFlag<std::string> label_list_path(parser, "LABEL_NUMBERS", "Specify labels to use in text file (one per line, default is all)", {'l', "labels"});
    args::ValueFlag<std::string> label_names_path(parser, "LABEL_NAMES", "Specify label names and print in first column (LABEL_NUMBERS must be specified)", {'n', "names"});
    args::Flag     transpose(parser, "TRANSPOSE", "Transpose output table (values go in rows instead of columns", {'t', "transpose"});
    args::Flag     ignore_zero(parser, "IGNORE_ZERO", "Ignore 0 label (background)", {'z', "ignore_zero"});
    args::ValueFlagList<std::string> header_paths(parser, "HEADER", "Add a header (can be specified multiple times)", {'H', "header"});
    QI::ParseArgs(parser, argc, argv);

    TLblGeoFilter::Pointer label_filter = TLblGeoFilter::New();
    label_filter->CalculatePixelIndicesOff();
    label_filter->CalculateOrientedBoundingBoxOff();
    label_filter->CalculateOrientedLabelRegionsOff();
    if (verbose) std::cout << "Reading first file: " << QI::CheckList(in_paths).at(0) << std::endl;
    QI::VolumeI::Pointer img = QI::ReadImage<QI::VolumeI>(QI::CheckList(in_paths).at(0));
    label_filter->SetInput(img);
    label_filter->Update();

    // Setup label number list
    typename TLblGeoFilter::LabelsType labels;
    std::vector<std::string> label_names;
    if (label_list_path) {
        if (verbose) std::cout << "Opening label numbers: " << label_list_path.Get() << std::endl;
        std::ifstream file(label_list_path.Get());
        std::string line;
        while (std::getline(file, line)) {
            std::stringstream linestream(line);
            int label;
            linestream >> label;
            labels.push_back(label);
        }
        if (label_names_path) {
            if (verbose) std::cout << "Opening label names: " << label_names_path.Get() << std::endl;
            std::ifstream name_file(label_names_path.Get());
            while (std::getline(name_file, line)) {
                label_names.push_back(line);
            }
            if (label_names.size() != labels.size()) {
                QI_EXCEPTION("Number of labels and names does not match");
            }
        }
    } else {
        if (verbose) std::cout << "Using default labels." << std::endl;
        labels = label_filter->GetLabels();
        std::sort(labels.begin(), labels.end());
    }
    if (ignore_zero) {
        if (verbose) std::cout << "Removing zero from label list." << std::endl;
        labels.erase(std::remove(labels.begin(), labels.end(), 0), labels.end());
    }

    // Setup headers (if any)
    std::vector<std::vector<std::string>> headers(header_paths.Get().size(), std::vector<std::string>());
    for (int h = 0; h < headers.size(); h++) {
        const std::string header_path = header_paths.Get().at(h);
        if (verbose) std::cout << "Reading header file: " << header_path << std::endl;
        std::ifstream header_file(header_path);
        std::string element;
        while (std::getline(header_file, element)) {
            headers.at(h).push_back(element);
        }
        std::cout << "h size: " << headers.at(h).size() << " files: " << in_paths.Get().size() << std::endl;
        if (headers.at(h).size() != in_paths.Get().size()) {
            QI_EXCEPTION("Number of input files (" << in_paths.Get().size() 
                         << ") does not match number of entries (" << headers.at(h).size()
                         << ") in header: " << header_path)
        }

    }

    double vox_volume = VoxelVolume(img);
    std::vector<std::vector<double>> volumes(in_paths.Get().size(), std::vector<double>(labels.size()));

    for (int i = 0; i < labels.size(); i++) {
        volumes.at(0).at(i) = vox_volume * label_filter->GetVolume(labels.at(i));
    }
    for (int f = 1; f < in_paths.Get().size(); f++) {
        if (verbose) std::cout << "Reading file: " << in_paths.Get().at(f) << std::endl;
        img = QI::ReadImage<QI::VolumeI>(in_paths.Get().at(f));
        label_filter->SetInput(img);
        label_filter->Update();
        vox_volume = VoxelVolume(img);
        for (int i = 0; i < labels.size(); i++) {
            volumes.at(f).at(i) = vox_volume * label_filter->GetVolume(labels.at(i));
        }
    }

    if (verbose) std::cout << "Writing CSV: " << std::endl;
    if (transpose) {
        if (label_names_path) {
            for (int i = 0; i < headers.size(); ++i) {
                std::cout << ",";
            }
            auto it = label_names.begin();
            std::cout << *it;
            for (++it; it != label_names.end(); ++it) {
                std::cout << "," << *it;
            }
            std::cout << std::endl;
        }
        for (int f = 0; f < volumes.size(); ++f){
            for (auto h = headers.begin(); h != headers.end(); h++) {
                std::cout << h->at(f) << ",";
            }
            auto v = volumes.at(f).begin();
            std::cout << *v;
            for (++v; v != volumes.at(f).end(); ++v) {
                std::cout << "," << *v;
            }
            std::cout << std::endl;
        }
    } else {
        for (auto hdr = headers.begin(); hdr != headers.end(); ++hdr) {
            if (label_names_path) std::cout << ",";
            auto h_el = hdr->begin();
            std::cout << *h_el;
            for (; h_el != hdr->end(); ++h_el) {
                std::cout << "," << *h_el;
            }
            std::cout << std::endl;
        }
        for (int l = 0; l < labels.size(); ++l) {
            if (label_names_path)
                std::cout << label_names.at(l) << ",";
            
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

