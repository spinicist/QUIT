/*
 *  qiaffine.cpp
 *
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>

#include "itkMetaDataObject.h"

#include "QI/Util.h"
#include "QI/IO.h"
#include "QI/Args.h"

using namespace std;

template<typename TImage>
void Convert(const std::string &input, const std::string &output) {
    auto image = QI::ReadImage<TImage>(input);
    QI::WriteImage<TImage>(image, output);
}

int main(int argc, char **argv) {
    args::ArgumentParser parser("Converts images between formats\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> input_file(parser, "INPUT", "Input file.");
    args::Positional<std::string> output_file(parser, "OUTPUT", "Output file.");

    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<std::string> out_prefix(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlagList<std::string> rename(parser, "RENAME", "Rename using specified header fields", {'r', "rename"});
    QI::ParseArgs(parser, argc, argv);

    const std::string input = QI::CheckPos(input_file);
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(input.c_str(), itk::ImageIOFactory::ReadMode);
    if (!imageIO) {
        cerr << "Could not open: " << input << endl;
        return EXIT_FAILURE;
    }
    imageIO->SetFileName(input);
    imageIO->ReadImageInformation();
    size_t dims = imageIO->GetNumberOfDimensions();
    auto PixelType = imageIO->GetPixelType();

    std::string output = out_prefix.Get();
    if (rename) {
        bool append_delim = false;
        for (const auto rename_field: args::get(rename)) {
            std::vector<std::string> string_array_value;
            std::string string_value;
            double double_value;
            auto dict = imageIO->GetMetaDataDictionary();
            if (!dict.HasKey(rename_field)) {
                std::cout << "Rename field '" << rename_field << "' not found in header. Ignoring" << std::endl;
                continue;
            }
            if (append_delim) {
                output.append("_");
            } else {
                append_delim = true;
            }
            if (ExposeMetaData(dict, rename_field, string_array_value)) {
                output.append(string_array_value[0]);
            } else if (ExposeMetaData(dict, rename_field, string_value)) {
                output.append(string_value);
            } else if (ExposeMetaData(dict, rename_field, double_value)) {
                std::ostringstream formatted;
                formatted << double_value;
                output.append(formatted.str());
            } else {
                QI_EXCEPTION("Could not determine type of rename header field:" << rename_field);
            }
        }
    } else {
        output.append(QI::Basename(input));
    }
    output.append(QI::CheckPos(output_file));
    #define DIM_SWITCH( N ) \
    case N:\
        switch (PixelType) {\
        case itk::ImageIOBase::SCALAR: Convert<itk::Image<float, N>>(input, output); break;\
        case itk::ImageIOBase::COMPLEX: Convert<itk::Image<std::complex<float>, N>>(input, output); break;\
        default: std::cerr << "Unsupported PixelType: " << PixelType << std::endl; return EXIT_FAILURE;\
        }\
    break;

    switch (dims) {
    DIM_SWITCH( 2 )
    DIM_SWITCH( 3 )
    DIM_SWITCH( 4 )
    }
    
    #undef DIM_SWITCH

    return EXIT_SUCCESS;
}
