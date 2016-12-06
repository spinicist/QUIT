/*
 *  qi_glmsetup.cpp
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
#include <algorithm>

#include "itkTileImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkDivideImageFilter.h"
#include "Filters/MeanImageFilter.h"

#include "QI/Types.h"
#include "QI/Util.h"
#include "QI/Option.h"

const std::string usage{
"Usage is: qi_glmsetup [options] files\n\
\n\
A utility for setting up merged 4D files for use with FSL randomise etc.\n\
Input individual files as arguments.\n\
The group for each input file must be specified with --groups (one per line)\n\
A value of 0 means ignore this file.\n\
An output file must be specified with --out.\n\
\n\
Additionally, simple GLM including co-variates can be generated.\n"};

int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::OptionList opts(usage);
    QI::Switch sort('s',"sort","Sort merged file and GLM in ascending group order", opts);
    QI::Option<std::string> covars_path("",'C',"covars","Path to covariates file (added to design matrix)", opts);
    QI::Option<std::string> ftests_path("",'f',"ftests","Generate and save F-tests", opts);
    QI::Option<std::string> contrasts_path("",'c',"contrasts","Generate and save contrasts", opts);
    QI::Option<std::string> design_path("",'d',"design","Path to save design matrix", opts);
    QI::Option<std::string> output_path("",'o',"out","Path for output merged file", opts);
    QI::Option<std::string> group_path("",'g',"groups","File to read group numbers from", opts);
    QI::Switch verbose('v',"verbose","Print more information", opts);
    QI::Help help(opts);
    std::deque<std::string> file_paths = opts.parse(argc, argv);
    if (!group_path.set()) {
        std::cerr << opts << std::endl;
        std::cerr << "Group file must be set with --groups option" << std::endl;
        return EXIT_FAILURE;
    }
    if (!output_path.set()) {
        std::cerr << opts << std::endl;
        std::cerr << "Output file must be set with --out option" << std::endl;
        return EXIT_FAILURE;
    }

    std::ifstream group_file(*group_path);
    if (!group_file) {
        std::cerr << "Group file: " << *group_path << " does not exist" << std::endl;
        return EXIT_FAILURE;
    }
    if (*verbose) std::cout << "Reading group file" << std::endl;
    std::vector<int> group_list;
    int temp;
    while (group_file >> temp) {
        group_list.push_back(temp);
    }
    int n_groups = *std::max_element(group_list.begin(), group_list.end());
    int n_images = std::count_if(group_list.begin(), group_list.end(), [](int i){return i > 0;}); // Count non-zero elements

    if (file_paths.size() != group_list.size()) {
        std::cerr << "Group list size and number of files do not match." << std::endl;
        return EXIT_FAILURE;
    }

    std::vector<std::vector<QI::VolumeF::Pointer>> groups(n_groups);
    if (*verbose) std::cout << "Number of groups: " << n_groups << std::endl;
    if (*verbose) std::cout << "Number of images: " << n_images << std::endl;
    itk::FixedArray<unsigned int, 4> layout;
    layout[0] = layout[1] = layout[2] = 1;
    layout[3] = n_images;
    auto tiler = itk::TileImageFilter<QI::VolumeF, QI::SeriesF>::New();
    tiler->SetLayout(layout);

    std::ofstream design_file;
    if (design_path.set()) {
        if (*verbose) std::cout << "Design matrix will be saved to: " << *design_path << std::endl;
        design_file = std::ofstream(*design_path);
    }
    std::vector<std::vector<std::vector<double>>> covars(n_groups);
    std::ifstream covars_file;
    if (covars_path.set()) {
        covars_file = std::ifstream(*covars_path);
    }
    int out_index = 0;
    for (int i = 0; i < group_list.size(); i++) {
        const int group = group_list.at(i);
        if (group > 0) { // Ignore entries with a 0
            if (*verbose) std::cout << "File: " << file_paths.at(i) << " Group: " << group << std::flush;
            QI::VolumeF::Pointer ptr = QI::ReadImage(file_paths.at(i));
            groups.at(group - 1).push_back(ptr);
            std::vector<double> covar;
            if (covars_path.set()) {
                std::string covar_line;
                std::getline(covars_file, covar_line);
                std::stringstream css(covar_line);
                double c;
                while (css >> c) {
                    covar.push_back(c);
                } 
                covars.at(group - 1).push_back(covar);
                if (*verbose) std::cout << ", read covariates.";
            }
            if (*verbose) std::cout << std::endl;
            if (!*sort) {
                tiler->SetInput(out_index, ptr);
                out_index++;
                if (design_path.set()) {
                    for (int g = 1; g <= n_groups; g++) {
                        if (g == group) {
                            design_file << "1\t";
                        } else {
                            design_file << "0\t";
                        }
                    }
                    if (covars_path.set()) {
                        for (auto& c : covar) {
                            design_file << c << "\t";
                        }
                    }
                    design_file << std::endl;
                }
            }
        } else {
            if (*verbose) std::cout << "Ignoring file: " << file_paths.at(i) << std::endl;
            // Eat a line in the covariates file as well
            if (covars_path.set()) {
                std::string covar_line;
                std::getline(covars_file, covar_line);
            }
        }
    }
    if (*sort) {
        if (*verbose) std::cout << "Sorting." << std::endl;
        for (int g = 0; g < n_groups; g++) {
            for (int i = 0; i < groups.at(g).size(); i++) {
                tiler->SetInput(out_index, groups.at(g).at(i));
                out_index++;
                if (design_path.set()) {
                    for (int g2 = 0; g2 < n_groups; g2++) {
                        if (g2 == g) {
                            design_file << "1\t";
                        } else {
                            design_file << "0\t";
                        }
                    }
                    if (covars_path.set()) {
                        for (auto& c : covars.at(g).at(i)) {
                            design_file << c << "\t";
                        }
                    }
                    design_file << std::endl;
                }
            }
        }
    }
    int n_covars = covars_path.set() ? covars.front().front().size() : 0;
    if (contrasts_path.set()) {
        if (*verbose) std::cout << "Generating contrasts" << std::endl;
        std::ofstream con_file(*contrasts_path);
        for (int g = 0; g < n_groups; g++) {
            for (int g2 = 0; g2 < n_groups; g2++) {
                if (g2 == g) {
                    con_file << "1\t";
                } else if (g2 == ((g+1) % (n_groups))) {
                    con_file << "-1\t";
                } else {
                    con_file << "0\t";
                }
            }
            for (int c = 0; c < n_covars; c++) {
                con_file << "0\t";
            }
            con_file << std::endl;
        }
        for (int c = 0; c < 2*n_covars; c++) {
            for (int g = 0; g < n_groups; g++) {
                con_file << "0\t";
            }
            for (int c2 = 0; c2 < n_covars; c2++) {
                if ((c/2) == c2) {
                    con_file << ((c % 2 == 0) ? "1\t" : "-1\t");
                } else {
                    con_file << "0\t";
                }
            }
            con_file << std::endl;
        }
    }
    if (ftests_path.set()) {
        if (*verbose) std::cout << "Generating F-tests" << std::endl;
        std::ofstream fts_file(*ftests_path);
        for (int g = 0; g < n_groups; g++) { // Individual group comparisons
            for (int g2 = 0; g2 < n_groups; g2++) {
                fts_file << ((g2 == g) ? "1\t" : "0\t");
            }
            for (int c = 0; c < 2*n_covars; c++) {
                fts_file << "0\t";
            }
            fts_file << std::endl;
        }
        if (n_groups > 2) { // Main effect
            for (int g = 0; g < (n_groups - 1); g++) {
                fts_file << "1\t";
            }
            fts_file << "0\t";
            for (int c = 0; c < 2*n_covars; c++) {
                fts_file << "0\t";
            }
            fts_file << std::endl;
        }
    }
    if (*verbose) std::cout << "Writing merged file: " << *output_path << std::endl;
    tiler->UpdateLargestPossibleRegion();
    // Reset space information because tiler messes it up
    QI::SeriesF::Pointer output = tiler->GetOutput();
    output->DisconnectPipeline();
    auto spacing    = output->GetSpacing();
    auto origin     = output->GetOrigin();
    auto direction  = output->GetDirection();
    spacing.Fill(1);
    origin.Fill(0);
    direction.SetIdentity();
    for (int i = 0; i < 3; i++) {
        spacing[i] = groups.at(0).at(0)->GetSpacing()[i];
        origin[i] =  groups.at(0).at(0)->GetOrigin()[i];
        for (int j = 0; j < 3; j++) {
            direction[i][j] = groups.at(0).at(0)->GetDirection()[i][j];
        }
    }
    output->SetSpacing(spacing);
    output->SetOrigin(origin);
    output->SetDirection(direction);
    QI::WriteImage<QI::SeriesF>(output, *output_path);
    return EXIT_SUCCESS;
}

