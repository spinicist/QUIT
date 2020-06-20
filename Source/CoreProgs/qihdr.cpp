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

#include "Args.h"
#include "ImageIO.h"
#include "Util.h"
#include "itkImageFileReader.h"
#include "itkMetaDataObject.h"

int hdr_main(args::Subparser &parser) {
    args::PositionalList<std::string> filenames(parser, "FILES", "Input files");

    args::Flag print_direction(
        parser, "DIRECTION", "Print the image direction/orientation", {'d', "direction"});
    args::Flag              print_origin(parser, "ORIGIN", "Print the the origin", {'o', "origin"});
    args::ValueFlag<size_t> print_spacing(
        parser, "SPACING", "Print voxel spacing (can specify one dimension)", {'S', "spacing"});
    args::ValueFlag<size_t> print_size(
        parser, "SIZE", "Print the matrix size (can specify one dimension)", {'s', "size"});
    args::Flag print_voxvol(
        parser, "VOLUME", "Calculate and print the volume of one voxel", {'v', "voxvol"});
    args::Flag print_type(parser, "DTYPE", "Print the data type", {'T', "dtype"});
    args::Flag print_dims(parser, "DIMS", "Print the number of dimensions", {'D', "dims"});
    args::Flag dim3(parser, "3D", "Treat input as 3D (discard higher dimensions)", {'3', "3D"});
    args::ValueFlagList<std::string> header_fields(
        parser,
        "METADATA",
        "Print a header metadata field (can be specified multiple times)",
        {'m', "meta"});
    parser.Parse();
    bool print_all = !(print_direction || print_origin || print_spacing || print_size ||
                       print_voxvol || print_type || print_dims || header_fields);
    for (const std::string &fname : QI::CheckList(filenames)) {
        itk::ImageIOBase::Pointer imageIO =
            itk::ImageIOFactory::CreateImageIO(fname.c_str(), itk::ImageIOFactory::ReadMode);
        if (!imageIO) {
            std::cerr << "Could not open: " << std::string(fname) << std::endl;
            break;
        }
        imageIO->SetFileName(std::string(fname));
        imageIO->ReadImageInformation();
        size_t dims = imageIO->GetNumberOfDimensions();
        if (print_all || verbose)
            fmt::print("File: {}\n", fname);
        if (print_all || verbose)
            fmt::print("Dimension:  ");
        if (print_all || print_dims)
            fmt::print("{}\n", dims);
        if (print_all || verbose)
            fmt::print("Voxel Type: ");
        if (print_all || print_type) {
            fmt::print("{} {}\n",
                       imageIO->GetPixelTypeAsString(imageIO->GetPixelType()),
                       imageIO->GetComponentTypeAsString(imageIO->GetComponentType()));
        }
        if (dim3 && dims > 3)
            dims = 3;
        if (print_all || print_size) {
            if (print_size.Get() > dims) {
                QI::Fail("Invalid dimension {} for image {}", print_size.Get(), fname);
            } else {
                int start_dim, end_dim;
                if (print_size.Get() == 0) {
                    start_dim = 0;
                    end_dim   = dims;
                } else {
                    start_dim = print_size.Get() - 1;
                    end_dim   = print_size.Get();
                }
                if (print_all || verbose)
                    std::cout << "Size:       ";
                if (print_all || print_size) {
                    for (int i = start_dim; i < end_dim; i++) {
                        std::cout << imageIO->GetDimensions(i);
                        if (i < (end_dim - 1))
                            std::cout << ",";
                    }
                    std::cout << std::endl;
                }
            }
        }
        if (print_all || print_spacing) {
            if (print_spacing.Get() > dims) {
                if (verbose)
                    std::cerr << "Invalid dimension " << print_spacing.Get() << " for image "
                              << fname << std::endl;
            } else {
                size_t start_dim, end_dim;
                if (print_spacing.Get() == 0) {
                    start_dim = 0;
                    end_dim   = dims;
                } else {
                    start_dim = print_spacing.Get() - 1;
                    end_dim   = print_spacing.Get();
                }
                if (print_all || verbose)
                    std::cout << "Spacing:    ";
                if (print_all || print_spacing) {
                    for (size_t i = start_dim; i < end_dim; i++)
                        std::cout << imageIO->GetSpacing(i) << "\t";
                    std::cout << std::endl;
                }
            }
        }
        if (print_all || verbose)
            std::cout << "Origin:     ";
        if (print_all || print_origin) {
            for (size_t i = 0; i < dims; i++)
                std::cout << imageIO->GetOrigin(i) << "\t";
            std::cout << std::endl;
        }
        if (print_all || verbose)
            std::cout << "Direction:  " << std::endl;
        if (print_all | print_direction) {
            for (size_t i = 0; i < dims; i++) {
                std::vector<double> dir = imageIO->GetDirection(i);
                for (size_t j = 0; j < dims; j++)
                    std::cout << dir[j] << "\t";
                std::cout << std::endl;
            }
        }
        if (print_all || verbose)
            std::cout << "Voxel vol:  ";
        if (print_all || print_voxvol) {
            double vol = imageIO->GetSpacing(0);
            for (size_t i = 1; i < dims; i++)
                vol *= imageIO->GetSpacing(i);
            std::cout << vol << std::endl;
        }
        for (const std::string &hf : header_fields.Get()) {
            auto header = imageIO->GetMetaDataDictionary();
            if (header.HasKey(hf)) {
                std::vector<std::string>              string_array_value;
                std::vector<std::vector<std::string>> string_array_array_value;
                std::vector<std::vector<double>>      double_array_array_value;
                std::vector<double>                   double_array_value;
                std::string                           string_value;
                double                                double_value;
                if (verbose)
                    std::cout << hf << ": ";
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
                    QI::Fail("Could not determine type of rename header field: {}", hf);
                }
            } else {
                QI::Log(verbose, "Header field not found: {}", hf);
            }
        }
    }
    return EXIT_SUCCESS;
}
