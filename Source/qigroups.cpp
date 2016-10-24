/*
 *  qigroups.cpp
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

#include "QI/Types.h"
#include "QI/Util.h"
#include "QI/Option.h"

const std::string usage{
"Usage is: qigroups [options] output_file\n\
\n\
A utility for setting up merged 4D files for use with FSL randomise etc.\n\
The list of files to be merged is read from stdin, one per line.\n\
A list of groups must be passed in with the --groups option. This\n\
file should contain a single number per line, indicating the group each\n\
matching line of stdin corresponds to. Group 0 is ignore.\n"
};

int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::OptionList opts(usage);
    QI::Option<std::string> group_path("",'g',"groups","File to read group numbers from", opts);
    QI::Option<std::string> output_path("",'o',"out","Path for output merged file", opts);
    QI::Option<std::string> design_path("",'d',"design","Path to save design matrix", opts);

    QI::Switch verbose('v',"verbose","Print more information", opts);
    QI::Help help(opts);
    std::vector<std::string> file_paths = opts.parse(argc, argv);
    if (!group_path.set()) {
        std::cerr << opts << std::endl;
        std::cerr << "Group file must be set with --groups option" << std::endl;
        return EXIT_FAILURE;
    }

    if (*verbose) std::cout << "Reading group file" << std::endl;
    std::ifstream group_file(*group_path);
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
    int out_index = 0;
    for (int i = 0; i < group_list.size(); i++) {
        const int group = group_list.at(i);
        if (group > 0) { // Ignore entries with a 0
            if (*verbose) std::cout << "Reading file: " << file_paths.at(i) << " for group " << group << std::endl;
            QI::VolumeF::Pointer ptr = QI::ReadImage(file_paths.at(i));
            groups.at(group - 1).push_back(ptr);
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
                design_file << std::endl;
            }
        } else {
            if (*verbose) std::cout << "Ignoring file: " << file_paths.at(i) << std::endl;
        }
    }
    if (*verbose) std::cout << "Writing merged file: " << *output_path << std::endl;
    tiler->UpdateLargestPossibleRegion();
    QI::WriteImage<QI::SeriesF>(tiler->GetOutput(), *output_path);
    return EXIT_SUCCESS;
}

