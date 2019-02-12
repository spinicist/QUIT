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

#include <algorithm>
#include <fstream>
#include <string>

#include "MeanImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkTileImageFilter.h"

#include "Args.h"
#include "ImageIO.h"
#include "ImageTypes.h"
#include "Util.h"

const std::string usage{"Usage is: qi_glmsetup [options] files\n\
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
    args::ArgumentParser parser(
        "A utility for setting up merged 4D files for use with FSL randomise etc.\n"
        "The group for each input file must be specified with --groups (one per line)\n"
        "A value of 0 means ignore this file.\n"
        "\nhttp://github.com/spinicist/QUIT");
    args::PositionalList<std::string> file_paths(parser, "INPUTS", "Input files to be merged.");
    args::HelpFlag                    help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag sort(parser, "SORT", "Sort merged file and GLM in ascending group order",
                    {'s', "sort"});
    args::ValueFlag<std::string> group_path(
        parser, "GROUPS", "File to read group numbers from (REQUIRED)", {'g', "groups"});
    args::ValueFlag<std::string> output_path(parser, "OUT", "Path for output merged file",
                                             {'o', "out"});
    args::ValueFlag<std::string> covars_path(
        parser, "COVARS", "Path to covariates file (added to design matrix)", {'C', "covars"});
    args::ValueFlag<std::string> design_path(parser, "DESIGN", "Path to save design matrix",
                                             {'d', "design"});
    args::ValueFlag<std::string> contrasts_path(parser, "CONTRASTS", "Generate and save contrasts",
                                                {'c', "contrasts"});
    args::ValueFlag<std::string> ftests_path(parser, "FTESTS", "Generate and save F-tests",
                                             {'f', "ftests"});
    QI::ParseArgs(parser, argc, argv, verbose);

    std::ifstream group_file(QI::CheckValue(group_path));
    QI::Log(verbose, "Reading group file");
    std::vector<int> group_list;
    int              temp;
    while (group_file >> temp) {
        group_list.push_back(temp);
    }
    int n_groups = *std::max_element(group_list.begin(), group_list.end());
    int n_images = std::count_if(group_list.begin(), group_list.end(),
                                 [](int i) { return i > 0; }); // Count non-zero elements

    if (file_paths.Get().size() != group_list.size()) {
        QI::Fail("Group list size and number of files do not match.");
    }

    std::vector<std::vector<QI::VolumeF::Pointer>> groups(n_groups);
    QI::Log(verbose, "Groups = {}, Images = {}", n_groups, n_images);
    itk::FixedArray<unsigned int, 4> layout;
    layout[0] = layout[1] = layout[2] = 1;
    layout[3]                         = n_images;
    auto tiler                        = itk::TileImageFilter<QI::VolumeF, QI::SeriesF>::New();
    tiler->SetLayout(layout);

    std::ofstream design_file;
    if (design_path) {
        QI::Log(verbose, "Design matrix will be saved to: {}", design_path.Get());
        design_file = std::ofstream(design_path.Get());
    }
    std::vector<std::vector<std::vector<std::string>>> covars(n_groups);
    std::vector<std::ifstream>                         covars_files;
    if (covars_path) {
        QI::Log(verbose, "All covariates: {}", covars_path.Get());
        std::istringstream stream_covars(covars_path.Get());
        while (!stream_covars.eof()) {
            std::string path;
            getline(stream_covars, path, ',');
            QI::Log(verbose, "Opening covariate file: {}", path);
            std::ifstream covars_file(path);
            if (!covars_file) {
                QI::Fail("Failed to open covariate file: {}", path);
            }
            covars_files.push_back(std::move(covars_file));
        }
    }
    int out_index = 0;
    for (size_t i = 0; i < group_list.size(); i++) {
        const int group = group_list.at(i);
        if (group > 0) { // Ignore entries with a 0
            QI::Log(verbose, "File: {} Group: {}", file_paths.Get().at(i), group);
            QI::VolumeF::Pointer ptr = QI::ReadImage(file_paths.Get().at(i), verbose);
            groups.at(group - 1).push_back(ptr);
            std::vector<std::string> covar;
            if (covars_path) {
                for (auto &f : covars_files) {
                    std::string c;
                    std::getline(f, c);
                    covar.push_back(c);
                }
                covars.at(group - 1).push_back(covar);
                QI::Log(verbose, "Covariate: {}", fmt::join(covar, ","));
            }
            if (!sort) {
                tiler->SetInput(out_index, ptr);
                out_index++;
                if (design_path) {
                    for (int g = 1; g <= n_groups; g++) {
                        if (g == group) {
                            design_file << "1\t";
                        } else {
                            design_file << "0\t";
                        }
                    }
                    if (covars_path) {
                        for (auto &c : covar) {
                            design_file << c << "\t";
                        }
                    }
                    design_file << "\n";
                }
            }
        } else {
            QI::Log(verbose, "Ignoring file: {}", file_paths.Get().at(i));
            // Eat a line in each covariates file
            for (auto &f : covars_files) {
                std::string dummy;
                std::getline(f, dummy);
            }
        }
    }
    if (sort) {
        QI::Log(verbose, "Sorting.");
        for (int g = 0; g < n_groups; g++) {
            for (size_t i = 0; i < groups.at(g).size(); i++) {
                tiler->SetInput(out_index, groups.at(g).at(i));
                out_index++;
                if (design_path) {
                    for (int g2 = 0; g2 < n_groups; g2++) {
                        if (g2 == g) {
                            design_file << "1\t";
                        } else {
                            design_file << "0\t";
                        }
                    }
                    if (covars_path) {
                        for (auto &c : covars.at(g).at(i)) {
                            design_file << c << "\t";
                        }
                    }
                    design_file << "\n";
                }
            }
        }
    }
    int n_covars = covars_path ? covars.front().front().size() : 0;
    if (contrasts_path) {
        QI::Log(verbose, "Generating contrasts");
        std::ofstream con_file(contrasts_path.Get());
        for (int g = 0; g < n_groups; g++) {
            for (int g2 = 0; g2 < n_groups; g2++) {
                if (g2 == g) {
                    con_file << "1\t";
                } else if (g2 == ((g + 1) % (n_groups))) {
                    con_file << "-1\t";
                } else {
                    con_file << "0\t";
                }
            }
            for (int c = 0; c < n_covars; c++) {
                con_file << "0\t";
            }
            con_file << "\n";
        }
        for (int c = 0; c < 2 * n_covars; c++) {
            for (int g = 0; g < n_groups; g++) {
                con_file << "0\t";
            }
            for (int c2 = 0; c2 < n_covars; c2++) {
                if ((c / 2) == c2) {
                    con_file << ((c % 2 == 0) ? "1\t" : "-1\t");
                } else {
                    con_file << "0\t";
                }
            }
            con_file << "\n";
        }
    }
    if (ftests_path) {
        QI::Log(verbose, "Generating F-tests");
        std::ofstream fts_file(ftests_path.Get());
        for (int g = 0; g < n_groups; g++) { // Individual group comparisons
            for (int g2 = 0; g2 < n_groups; g2++) {
                fts_file << ((g2 == g) ? "1\t" : "0\t");
            }
            for (int c = 0; c < 2 * n_covars; c++) {
                fts_file << "0\t";
            }
            fts_file << std::endl;
        }
        if (n_groups > 2) { // Main effect
            for (int g = 0; g < (n_groups - 1); g++) {
                fts_file << "1\t";
            }
            fts_file << "0\t";
            for (int c = 0; c < 2 * n_covars; c++) {
                fts_file << "0\t";
            }
            fts_file << std::endl;
        }
    }
    tiler->UpdateLargestPossibleRegion();
    // Reset space information because tiler messes it up
    QI::SeriesF::Pointer output = tiler->GetOutput();
    output->DisconnectPipeline();
    auto spacing   = output->GetSpacing();
    auto origin    = output->GetOrigin();
    auto direction = output->GetDirection();
    spacing.Fill(1);
    origin.Fill(0);
    direction.SetIdentity();
    for (int i = 0; i < 3; i++) {
        spacing[i] = groups.at(0).at(0)->GetSpacing()[i];
        origin[i]  = groups.at(0).at(0)->GetOrigin()[i];
        for (int j = 0; j < 3; j++) {
            direction[i][j] = groups.at(0).at(0)->GetDirection()[i][j];
        }
    }
    output->SetSpacing(spacing);
    output->SetOrigin(origin);
    output->SetDirection(direction);
    QI::WriteImage<QI::SeriesF>(output, QI::CheckValue(output_path), verbose);
    return EXIT_SUCCESS;
}
