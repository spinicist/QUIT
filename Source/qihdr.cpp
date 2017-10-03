/*
 *  qihdr.cpp
 *
 *  Created by Tobias Wood on 11/06/2015.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include "itkMetaDataObject.h"
#include "QI/IO.h"
#include "QI/Args.h"
#include "QI/Util.h"

args::ArgumentParser parser(
"Extracts information from image headers.\n"
"By default, a summary of the header is printed. If any options are specified,"
"only those parts of the header will be printed. Multiple files can be input"
"in which case the header info is written for each in order.\n"
"http://github.com/spinicist/QUIT");

args::PositionalList<std::string> filenames(parser, "FILES", "Input files");
args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
args::Flag     verbose(parser, "VERBOSE", "Print description of each output line", {'v', "verbose"});
args::Flag print_direction(parser, "DIRECTION", "Print the image direction/orientation", {'d', "direction"});
args::Flag print_origin(parser, "ORIGIN", "Print the the origin", {'o', "origin"});
args::ValueFlag<int> print_spacing(parser, "SPACING", "Print voxel spacing (can specify one dimension)", {'S',"spacing"});
args::ValueFlag<int> print_size(parser, "SIZE", "Print the matrix size (can specify one dimension)", {'s',"size"});
args::Flag print_voxvol(parser, "VOLUME", "Calculate and print the volume of one voxel", {'v',"voxvol"});
args::Flag print_type(parser, "DTYPE", "Print the data type",{'T',"dtype"});
args::Flag print_dims(parser, "DIMS", "Print the number of dimensions",{'D',"dims"});
args::Flag dim3(parser, "3D", "Treat input as 3D (discard higher dimensions)", {'3',"3D"});
args::ValueFlagList<std::string> header_fields(parser, "METADATA", "Print a header metadata field (can be specified multiple times)", {'m', "meta"});

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {

    QI::ParseArgs(parser, argc, argv);
    bool print_all = !(print_direction || print_origin || print_spacing || print_size ||
                       print_voxvol || print_type || print_dims || header_fields);
    for (const std::string& fname : QI::CheckList(filenames)) {
        itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(fname.c_str(), itk::ImageIOFactory::ReadMode);
        if (!imageIO) {
            std::cerr << "Could not open: " << std::string(fname) << std::endl;
            break;
        }
        imageIO->SetFileName(std::string(fname));
        imageIO->ReadImageInformation();
        size_t dims = imageIO->GetNumberOfDimensions();
        if (verbose) std::cout << "File:       " << std::string(fname) << std::endl;
        if (print_all || verbose) std::cout << "Dimension:  "; if (print_all || print_dims) std::cout << dims << std::endl;
        if (print_all || verbose) std::cout << "Voxel Type: ";
        if (print_all || print_type) {
            std::cout << imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) << " " 
                      << imageIO->GetComponentTypeAsString(imageIO->GetComponentType()) << std::endl;
        }
        if (dim3 && dims > 3) dims = 3;
        if (print_all || print_size) {
            if (print_size.Get() < 0 || print_size.Get() > dims) {
                if (verbose) std::cerr << "Invalid dimension " << print_size.Get() << " for image " << fname << std::endl;
            } else {
                int start_dim, end_dim;
                if (print_size.Get() == 0) {
                    start_dim = 0; end_dim = dims;
                } else {
                    start_dim = print_size.Get() - 1;
                    end_dim = print_size.Get();
                }
                if (print_all || verbose) std::cout << "Size:       "; if (print_all || print_size) { for (int i = start_dim; i < end_dim; i++) std::cout << imageIO->GetDimensions(i) << "\t"; std::cout << std::endl; }
            }
        }
        if (print_all || print_spacing) {
            if (print_spacing.Get() < 0 || print_spacing.Get() > dims) {
                if (verbose) std::cerr << "Invalid dimension " << print_spacing.Get() << " for image " << fname << std::endl;
            } else {
                int start_dim, end_dim;
                if (print_spacing.Get() == 0) {
                    start_dim = 0; end_dim = dims;
                } else {
                    start_dim = print_spacing.Get() - 1;
                    end_dim = print_spacing.Get();
                }
                if (print_all || verbose) std::cout << "Spacing:    "; if (print_all || print_spacing) { for (int i = start_dim; i < end_dim; i++) std::cout << imageIO->GetSpacing(i) << "\t"; std::cout << std::endl; }
            }
        }
        if (print_all || verbose) std::cout << "Origin:     "; if (print_all || print_origin)   { for (int i = 0; i < dims; i++) std::cout << imageIO->GetOrigin(i) << "\t"; std::cout << std::endl; }
        if (print_all || verbose) std::cout << "Direction:  " << std::endl;
        if (print_all | print_direction) {
            for (int i = 0; i < dims; i++) {
                std::vector<double> dir = imageIO->GetDirection(i);
                for (int j = 0; j < dims; j++)
                    std::cout << dir[j] << "\t";
                std::cout << std::endl;
            }
        }
        if (print_all || verbose) std::cout << "Voxel vol:  "; if (print_all || print_voxvol)   { double vol = imageIO->GetSpacing(0); for (int i = 1; i < dims; i++) vol *= imageIO->GetSpacing(i); std::cout << vol << std::endl; }
        for (const std::string &hf : header_fields.Get()) {
            auto header = imageIO->GetMetaDataDictionary();
            if (header.HasKey(hf)) {
                std::vector<std::string> string_array_value;
                std::vector<std::vector<std::string> > string_array_array_value;
                std::vector<std::vector<double> > double_array_array_value;
                std::vector<double> double_array_value;
                std::string string_value;
                auto delim = "";
                double double_value;
                if (verbose) std::cout << hf << ": ";
                if (ExposeMetaData(header, hf, string_array_value)) {
                    std::cout << string_array_value << std::endl;
                } else if (ExposeMetaData(header, hf, string_array_array_value)) {
                    std::cout << string_array_array_value << std::endl;
                } else if (ExposeMetaData(header, hf, double_array_value)) {
                    std::cout << double_array_value << std::endl;
                } else if (ExposeMetaData(header, hf, double_array_array_value)) {
                    std::cout << double_array_array_value << std::endl;
                } else if (ExposeMetaData(header, hf, string_value)) {
                   std::cout << string_value << std::endl;
                } else if (ExposeMetaData(header, hf, double_value)) {
                    std::cout << double_value << std::endl;
                } else {
                    QI_EXCEPTION("Could not determine type of rename header field:" << hf);
                }
            } else {
                if (verbose) std::cout << "Header field not found: " << hf << std::endl;
            }
        }
    }
    return EXIT_SUCCESS;
}
