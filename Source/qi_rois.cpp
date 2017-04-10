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

// Declare arguments here so they are available in helper functions
args::ArgumentParser parser("Calculates average values or volumes of ROI labels.\n"
                            "If the --volumes flag is specified, only give the label images.\n"
                            "If --volumes is not specified, give label and value image pairs (all labels first, then all values).\n"
                            "Output goes to stdout\n"
                            "http://github.com/spinicist/QUIT");
args::PositionalList<std::string> in_paths(parser, "INPUT", "Input file paths.");
args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
args::Flag     verbose(parser, "VERBOSE", "Print more information (will mess up output, for test runs only)", {'v', "verbose"});
args::Flag     volumes(parser, "VOLUMES", "Output ROI volumes, not average values (does not require value images)", {'V', "volumes"});
args::ValueFlag<std::string> label_list_path(parser, "LABEL_NUMBERS", "Specify labels to use in text file (one per line, default is all)", {'l', "labels"});
args::ValueFlag<std::string> label_names_path(parser, "LABEL_NAMES", "Specify label names and print in first column (LABEL_NUMBERS must be specified)", {'n', "names"});
args::Flag     transpose(parser, "TRANSPOSE", "Transpose output table (values go in rows instead of columns", {'t', "transpose"});
args::Flag     ignore_zero(parser, "IGNORE_ZERO", "Ignore 0 label (background)", {'z', "ignore_zero"});
args::ValueFlagList<std::string> header_paths(parser, "HEADER", "Add a header (can be specified multiple times)", {'H', "header"});

/*
 * Helper function to calculate the volume of a voxel in an image
 */
double VoxelVolume(QI::VolumeI::Pointer img) {
    double vox_volume = img->GetSpacing()[0];
    for (int v = 1; v < 3; v++)
        vox_volume *= img->GetSpacing()[v];
    return vox_volume;
}

/*
 * Helper function to work out the label list
 */
void GetLabelList(TLblGeoFilter::LabelsType &labels, std::vector<std::string> &label_names) {
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
        if (verbose) std::cout << "Reading first label file to determine labels: " << QI::CheckList(in_paths).at(0) << std::endl;
        TLblGeoFilter::Pointer label_filter = TLblGeoFilter::New();
        label_filter->CalculatePixelIndicesOff();
        label_filter->CalculateOrientedBoundingBoxOff();
        label_filter->CalculateOrientedLabelRegionsOff();
        QI::VolumeI::Pointer img = QI::ReadImage<QI::VolumeI>(QI::CheckList(in_paths).at(0));
        label_filter->SetInput(img);
        label_filter->Update();
        labels = label_filter->GetLabels();
        std::sort(labels.begin(), labels.end());
    }
    if (ignore_zero) {
        if (verbose) std::cout << "Removing zero from label list." << std::endl;
        labels.erase(std::remove(labels.begin(), labels.end(), 0), labels.end());
    }
}

/*
 * Helper function to read in all the header lines
 */
std::vector<std::vector<std::string>> GetHeaders(int n_files) {
    std::vector<std::vector<std::string>> headers(header_paths.Get().size(), std::vector<std::string>());
    for (int h = 0; h < headers.size(); h++) {
        const std::string header_path = header_paths.Get().at(h);
        if (verbose) std::cout << "Reading header file: " << header_path << std::endl;
        std::ifstream header_file(header_path);
        std::string element;
        while (std::getline(header_file, element)) {
            headers.at(h).push_back(element);
        }
        if (headers.at(h).size() != n_files) {
            QI_EXCEPTION("Number of input files (" << n_files 
                         << ") does not match number of entries (" << headers.at(h).size()
                         << ") in header: " << header_path)
        }
    }
    return headers;
}

/*
 * Helper function to actually work out all the values
 */
std::vector<std::vector<double>> GetValues(const int n_files, const TLblGeoFilter::LabelsType &labels) {
    std::vector<std::vector<double>> values(in_paths.Get().size(), std::vector<double>(labels.size()));
    TLblGeoFilter::Pointer label_filter = TLblGeoFilter::New();
    QI::VolumeI::Pointer label_img = ITK_NULLPTR;
    QI::VolumeF::Pointer value_img = ITK_NULLPTR;
    for (int f = 0; f < n_files; f++) {
        if (volumes) {
            if (verbose) std::cout << "Reading label file: " << in_paths.Get().at(f) << std::endl;
            label_img = QI::ReadImage<QI::VolumeI>(in_paths.Get().at(f));
        } else {
            if (verbose) std::cout << "Reading label file: " << in_paths.Get().at(f) << std::endl;
            label_img = QI::ReadImage<QI::VolumeI>(in_paths.Get().at(f));
            if (verbose) std::cout << "Reading value file: " << in_paths.Get().at(f + n_files) << std::endl;
            value_img = QI::ReadImage(in_paths.Get().at(f + n_files));
            label_filter->SetIntensityInput(value_img);
        }
        double vox_volume = VoxelVolume(label_img);
        label_filter->SetInput(label_img);
        label_filter->Update();
        for (int i = 0; i < labels.size(); i++) {
            const double volume = vox_volume * label_filter->GetVolume(labels.at(i));
            if (volumes) {
                values.at(f).at(i) = volume;
            } else {
                const double integrated = label_filter->GetIntegratedIntensity(labels.at(i));
                values.at(f).at(i) = integrated / volume;
            }
        }
    }
    return values;
}

/*
 * MAIN
 */
int main(int argc, char **argv) {
    QI::ParseArgs(parser, argc, argv);

    int n_files = QI::CheckList(in_paths).size();
    if (volumes) {
        if (verbose) std::cout << "There are " << n_files << " input images, finding ROI volumes" << std::endl;
    } else {
        // For ROI values, we have label and value image pairs
        if (n_files % 2 != 0) {
            std::cerr << "Require an even number of input images when finding ROI values" << std::endl;
            return EXIT_FAILURE;
        }
        n_files = n_files / 2;
        if (verbose) std::cout << "There are " << n_files << " input image pairs, finding mean ROI values" << std::endl;
    }
    // Setup label number list
    typename TLblGeoFilter::LabelsType labels;
    std::vector<std::string> label_names;
    GetLabelList(labels, label_names);

    // Setup headers (if any)
    auto headers = GetHeaders(n_files);

    // Now get the values/volumes
    auto values_table = GetValues(n_files, labels);

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
        for (int row = 0; row < values_table.size(); ++row){
            for (auto h = headers.begin(); h != headers.end(); h++) {
                std::cout << h->at(row) << ",";
            }
            auto values_row_it = values_table.at(row).begin();
            std::cout << *values_row_it;
            for (++values_row_it; values_row_it != values_table.at(row).end(); ++values_row_it) {
                std::cout << "," << *values_row_it;
            }
            std::cout << std::endl;
        }
    } else {
        for (auto hdr = headers.begin(); hdr != headers.end(); ++hdr) {
            if (label_names_path) std::cout << ",";
            auto h_el = hdr->begin();
            std::cout << *h_el;
            for (++h_el; h_el != hdr->end(); ++h_el) {
                std::cout << "," << *h_el;
            }
            std::cout << std::endl;
        }
        for (int l = 0; l < labels.size(); ++l) {
            if (label_names_path)
                std::cout << label_names.at(l) << ",";
            
            auto values_col_it = values_table.begin();
            std::cout << values_col_it->at(l);
            for (++values_col_it; values_col_it != values_table.end(); ++values_col_it) {
                std::cout << "," << values_col_it->at(l);
            }
            std::cout << std::endl;
        }
    }
    return EXIT_SUCCESS;
}

