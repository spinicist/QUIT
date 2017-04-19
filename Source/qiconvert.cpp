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

template<int D, typename T>
void ConvertPixel(const std::string &input, const std::string &output,
                  const itk::ImageIOBase::IOPixelType &pix)
{
    switch (pix) {
        case itk::ImageIOBase::SCALAR: Convert<itk::Image<T, D>>(input, output); break;
        case itk::ImageIOBase::COMPLEX: Convert<itk::Image<std::complex<T>, D>>(input, output); break;
        default: QI_FAIL("Unsupported pixel type in image " << input);
    }
}

template<int D>
void ConvertComponent(const std::string &input, const std::string &output,
                      const itk::ImageIOBase::IOComponentType &component,
                      const itk::ImageIOBase::IOPixelType &pix)
{
    switch (component) {
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE: QI_FAIL("Unknown component type in image " << input);
        case itk::ImageIOBase::UCHAR:  ConvertPixel<D, unsigned char>(input, output, pix); break;
        case itk::ImageIOBase::CHAR:   ConvertPixel<D, char>(input, output, pix); break;
        case itk::ImageIOBase::USHORT: ConvertPixel<D, unsigned short>(input, output, pix); break;
        case itk::ImageIOBase::SHORT:  ConvertPixel<D, short>(input, output, pix); break;
        case itk::ImageIOBase::UINT:   ConvertPixel<D, unsigned int>(input, output, pix); break;
        case itk::ImageIOBase::INT:    ConvertPixel<D, int>(input, output, pix); break;
        case itk::ImageIOBase::ULONG:  ConvertPixel<D, unsigned long>(input, output, pix); break;
        case itk::ImageIOBase::LONG:   ConvertPixel<D, long>(input, output, pix); break;
        case itk::ImageIOBase::FLOAT:  ConvertPixel<D, float>(input, output, pix); break;
        case itk::ImageIOBase::DOUBLE: ConvertPixel<D, double>(input, output, pix); break;
    }
}

void ConvertDims(const std::string &input, const std::string &output, const int D,
                 const itk::ImageIOBase::IOComponentType &component,
                 const itk::ImageIOBase::IOPixelType &pix)
{
    switch (D) {
        case 2: ConvertComponent<2>(input, output, component, pix); break;
        case 3: ConvertComponent<3>(input, output, component, pix); break;
        case 4: ConvertComponent<4>(input, output, component, pix); break;
        default: QI_FAIL("Unsupported dimension: " << D);
    }
}

int main(int argc, char **argv) {
    args::ArgumentParser parser("Converts images between formats\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> input_file(parser, "INPUT", "Input file.");
    args::Positional<std::string> output_suffix(parser, "OUTPUT", "Output suffix, must include extension.");

    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<std::string> output_prefix(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlagList<std::string> rename(parser, "RENAME", "Rename using specified header fields", {'r', "rename"});
    QI::ParseArgs(parser, argc, argv);

    const std::string input = QI::CheckPos(input_file);
    itk::ImageIOBase::Pointer header = itk::ImageIOFactory::CreateImageIO(input.c_str(), itk::ImageIOFactory::ReadMode);
    if (!header) QI_FAIL("Could not open: " << input);
    header->SetFileName(input);
    header->ReadImageInformation();

    /* Deal with renaming */
    std::string output = output_prefix.Get();
    if (rename) {
        bool append_delim = false;
        for (const auto rename_field: args::get(rename)) {
            std::vector<std::string> string_array_value;
            std::string string_value;
            double double_value;
            auto dict = header->GetMetaDataDictionary();
            if (!dict.HasKey(rename_field)) {
                if (verbose) std::cout << "Rename field '" << rename_field << "' not found in header. Ignoring" << std::endl;
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
    output.append(QI::CheckPos(output_suffix));
    std::cout << output << std::endl;

    auto dims = header->GetNumberOfDimensions();
    auto pixel_type = header->GetPixelType();
    auto component_type = header->GetComponentType();

    ConvertDims(input, output, dims, component_type, pixel_type);

    return EXIT_SUCCESS;
}
