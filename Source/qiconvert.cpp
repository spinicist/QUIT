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
    Eigen::initParallel();
    QI::ArgParser args{argc, argv,
        "Usage is: qiconvert input output_format [options]\n"
        "Converts input to specified format",
        {{"help", 'h', "Display the help message and quit", false},
         {"verbose", 'v', "Print more information", false},
         {"prefix", 'p', "Prefix output path", true},
         {"rename", 'r', "Rename using specified header field", true}}
    };

    bool verbose = args.option_present("verbose");
    std::deque<const std::string> nonopts = args.nonoptions();
    if ((nonopts.size() == 0) || (nonopts.size() > 2)) {
        std::cerr << "Incorrect number of arguments, use -h to see usage." << std::endl;
        return EXIT_FAILURE;
    }

    const std::string input = nonopts[0];
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(input.c_str(), itk::ImageIOFactory::ReadMode);
    if (!imageIO) {
        cerr << "Could not open: " << input << endl;
    }
    imageIO->SetFileName(input);
    imageIO->ReadImageInformation();
    size_t dims = imageIO->GetNumberOfDimensions();
    auto PixelType = imageIO->GetPixelType();

    std::string output = args.option_value<string>("prefix","");
    if (args.option_present("rename")) {
        std::string rename_field = args.option_value<string>("rename","");
        std::vector<std::string> string_array_value;
        std::string string_value;
        double double_value;

        auto dict = imageIO->GetMetaDataDictionary();
        if (ExposeMetaData(dict, rename_field, string_array_value)) {
            output += string_array_value[0] + nonopts[1];
        } else if (ExposeMetaData(dict, rename_field, string_value)) {
            output += string_value + nonopts[1];
        } else if (ExposeMetaData(dict, rename_field, double_value)) {
            output += std::to_string(double_value) + nonopts[1];
        } else {
            QI_EXCEPTION("Could not determine type of rename header field");
        }
    } else {
        output += QI::StripExt(input) + nonopts[1];
    }

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
