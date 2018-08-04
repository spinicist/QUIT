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
#include <iomanip>
#include <fstream>
#include <string>

#include "ImageTypes.h"
#include "ImageIO.h"
#include "Args.h"
#include "Util.h"

#include "itkLabelStatisticsImageFilter.h"
typedef itk::LabelStatisticsImageFilter<QI::VolumeF, QI::VolumeI> TStatsFilter;
typedef TStatsFilter::ValidLabelValuesContainerType TLabels;
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
args::ValueFlag<std::string> label_list_path(parser, "LABELS", "Specify labels and names to use in text file 'label number, label name' one per line", {'l', "labels"});
args::Flag     print_names(parser, "PRINT_NAMES", "Print label names in first column/row (LABEL_NUMBERS must be specified)", {'n', "print_names"});
args::Flag     transpose(parser, "TRANSPOSE", "Transpose output table (values go in rows instead of columns", {'t', "transpose"});
args::Flag     ignore_zero(parser, "IGNORE_ZERO", "Ignore 0 label (background)", {'z', "ignore_zero"});
args::Flag     sigma(parser, "SIGMA", "Print ±std along with mean", {'s', "sigma"});
args::ValueFlag<int> precision(parser, "PRECISION", "Number of decimal places (default 6)", {'p', "precision"}, 6);
args::ValueFlag<std::string> delim(parser, "DELIMITER", "Specify delimiter to use between entries (default ,)", {'d',"delim"}, ",");
args::ValueFlagList<std::string> header_paths(parser, "HEADER", "Add a header (can be specified multiple times)", {'H', "header"});
args::ValueFlagList<std::string> header_names(parser, "HEADER NAME", "Header name (must be specified in same order as paths)", {"header_name"});
args::ValueFlagList<double> scales(parser, "SCALE", "Divide ROI values by scale (must be same order as paths)", {"scale"});

/*
 * Helper function to work out the label list
 */
void GetLabelList(TLabels &label_numbers, std::vector<std::string> &label_names) {
    if (label_list_path) {
        QI_LOG(verbose, "Opening label list file: " << label_list_path.Get());
        std::ifstream file(label_list_path.Get());
        if (!file) {
            QI_EXCEPTION("Could not open label list file: " << label_list_path.Get());
        }
        std::string temp;
        while (std::getline(file, temp, ',')) {
            label_numbers.push_back(stoi(temp));
            std::getline(file, temp);
            temp.erase(std::remove(temp.begin(), temp.end(), '\r'), temp.end()); // Deal with rogue ^M characters
            label_names.push_back(temp);
            QI_LOG(verbose, "Read label: " << label_numbers.back() << ", name: " << label_names.back());
        }
    } else {
        TStatsFilter::Pointer label_filter = TStatsFilter::New();
        QI_LOG(verbose, "Reading first label file to determine labels: " << QI::CheckList(in_paths).at(0));
        label_filter->SetLabelInput(QI::ReadImage<QI::VolumeI>(QI::CheckList(in_paths).at(0)));
        label_filter->SetInput(QI::ReadImage<QI::VolumeF>(QI::CheckList(in_paths).at(0))); // Dummy value because ITK complains input primary is not set
        label_filter->Update();
        label_numbers = label_filter->GetValidLabelValues();
        std::sort(label_numbers.begin(), label_numbers.end());
        if (verbose) {
            std::cout << "Found the following labels:\n";
            for (auto &l : label_numbers) std::cout << l << " ";
            std::cout << std::endl;
        }
    }
    if (ignore_zero) {
        QI_LOG(verbose, "Removing zero from label list.");
        label_numbers.erase(std::remove(label_numbers.begin(), label_numbers.end(), 0), label_numbers.end());
    }
}

/*
 * Helper function to read in all the header lines
 */
std::vector<std::vector<std::string>> GetHeaders(size_t n_files) {
    if (header_names && (header_names.Get().size() != header_paths.Get().size())) {
        QI_EXCEPTION("Number of header names must match number of header paths");
    }
    std::vector<std::vector<std::string>> headers(header_paths.Get().size(), std::vector<std::string>());
    for (size_t h = 0; h < headers.size(); h++) {
        const std::string header_path = header_paths.Get().at(h);
        QI_LOG(verbose, "Reading header file: " << header_path);
        std::ifstream header_file(header_path);
        if (!header_file) {
            QI_EXCEPTION("Failed to open header file: " << header_path);
        }
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
void GetValues(const int n_files, const TLabels &labels, const std::vector<double> &scale_list,
               std::vector<std::vector<double>> &mean_table, std::vector<std::vector<double>> &sigma_table, std::vector<std::vector<double>> &volume_table) {
    mean_table = std::vector<std::vector<double>>(n_files, std::vector<double>(labels.size()));
    sigma_table = std::vector<std::vector<double>>(n_files, std::vector<double>(labels.size()));
    volume_table = std::vector<std::vector<double>>(n_files, std::vector<double>(labels.size()));
    TStatsFilter::Pointer label_filter = TStatsFilter::New();
    QI::VolumeI::Pointer label_img = ITK_NULLPTR;
    QI::VolumeF::Pointer value_img = ITK_NULLPTR;
    for (int f = 0; f < n_files; f++) {
        if (volumes) {
            QI_LOG(verbose, "Reading label file: " << in_paths.Get().at(f));
            label_img = QI::ReadImage<QI::VolumeI>(in_paths.Get().at(f));
            value_img = QI::ReadImage(in_paths.Get().at(f)); // Dummy image
        } else {
            QI_LOG(verbose, "Reading label file: " << in_paths.Get().at(f));
            label_img = QI::ReadImage<QI::VolumeI>(in_paths.Get().at(f));
            QI_LOG(verbose, "Reading value file: " << in_paths.Get().at(f + n_files));
            value_img = QI::ReadImage(in_paths.Get().at(f + n_files));
        }
        double vox_volume = QI::VoxelVolume(label_img);
        label_filter->SetLabelInput(label_img);
        label_filter->SetInput(value_img);
        label_filter->Update();
        for (size_t i = 0; i < labels.size(); i++) {
            mean_table.at(f).at(i) = label_filter->GetMean(labels.at(i)) / scale_list.at(f);
            sigma_table.at(f).at(i) = label_filter->GetSigma(labels.at(i)) / scale_list.at(f);
            volume_table.at(f).at(i) = label_filter->GetCount(labels.at(i)) * vox_volume;
        }
    }
}

/*
 * MAIN
 */
int main(int argc, char **argv) {
    QI::ParseArgs(parser, argc, argv, verbose);

    size_t n_files = QI::CheckList(in_paths).size();
    if (volumes) {
        QI_LOG(verbose, "There are " << n_files << " input images, finding ROI volumes");
    } else {
        // For ROI values, we have label and value image pairs
        if (n_files % 2 != 0) {
            QI_FAIL("Require an even number of input images when finding ROI values");
        }
        n_files = n_files / 2;
        QI_LOG(verbose, "There are " << n_files << " input image pairs, finding mean ROI values");
    }
    // Setup label number list
    TLabels labels;
    std::vector<std::string> label_names;
    GetLabelList(labels, label_names);

    // Setup headers (if any)
    auto headers = GetHeaders(n_files);

    // Setup scales
    std::vector<double> scale_list = scales.Get();
    if (scale_list.size() < n_files) {
        int old_size = scale_list.size();
        scale_list.resize(n_files);
        for (size_t i = old_size; i < n_files; i++) {
            scale_list.at(i) = 1.0;
        }
    }
    if (verbose) {
        std::cout << "Scales are: ";
        for (const auto &s: scale_list) std::cout << s << " ";
        std::cout << std::endl;
    }

    // Now get the values/volumes
    std::vector<std::vector<double>> mean_table, sigma_table, volume_table;
    GetValues(n_files, labels, scale_list, mean_table, sigma_table, volume_table);

    QI_LOG(verbose, "Writing CSV: ");
    if (precision) std::cout << std::fixed << std::setprecision(precision.Get());
    if (transpose) {
        if (print_names) {
            for (size_t h = 0; h < headers.size(); ++h) {
                if (header_names) {
                    std::cout << header_names.Get().at(h) << delim.Get();
                } else {
                    std::cout << header_paths.Get().at(h) << delim.Get();
                }
            }
            auto it = label_names.begin();
            std::cout << *it;
            for (++it; it != label_names.end(); ++it) {
                std::cout << delim.Get() << *it;
            }
            std::cout << std::endl;
        }
        for (size_t row = 0; row < n_files; ++row){
            for (auto h = headers.begin(); h != headers.end(); h++) {
                std::cout << h->at(row) << delim.Get();
            }
            if (volumes) {
                std::cout << volume_table.at(row).at(0);
            } else {
                std::cout << mean_table.at(row).at(0);
                if (sigma) std::cout << "±" << sigma_table.at(row).at(0);
            }
            for (size_t val = 1; val < labels.size(); ++val) {
                std::cout << delim.Get();
                if (volumes) {
                    std::cout << volume_table.at(row).at(val);
                } else {
                    std::cout << mean_table.at(row).at(val);
                    if (sigma) std::cout << "±" << sigma_table.at(row).at(val);
                }
            }
            std::cout << std::endl;
        }
    } else {
        for (size_t h = 0; h < headers.size(); h++) {
            if (print_names) {
                if (header_names) {
                    std::cout << header_names.Get().at(h) << delim.Get();
                } else {
                    std::cout << header_paths.Get().at(h) << delim.Get();
                }
            }
            auto h_el = headers.at(h).begin();
            std::cout << *h_el;
            for (++h_el; h_el != headers.at(h).end(); ++h_el) {
                std::cout << delim.Get() << *h_el;
            }
            std::cout << std::endl;
        }
        for (size_t l = 0; l < labels.size(); ++l) {
            if (print_names) {
                std::cout << label_names.at(l) << std::flush << delim.Get();
            }

            if (volumes) {
                std::cout << volume_table.at(0).at(l);
            } else {
                std::cout << mean_table.at(0).at(l);
                if (sigma) std::cout << "±" << sigma_table.at(0).at(l);
            }
            for (size_t file = 1; file < n_files; ++file) {
                std::cout << delim.Get();
                if (volumes) {
                    std::cout << volume_table.at(file).at(l);
                } else {
                    std::cout << mean_table.at(file).at(l);
                    if (sigma) std::cout << "±" << sigma_table.at(0).at(l);
                }
            }
            std::cout << std::endl;
        }
    }
    return EXIT_SUCCESS;
}

