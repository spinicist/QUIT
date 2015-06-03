/*
 *  niicomplex.cpp
 *
 *  Created by Tobias Wood on 03/06/2015.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <getopt.h>
#include <iostream>
#include <random>
#include <functional>

#include "itkImage.h"
#include "itkComplexToPhaseImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkComposeImageFilter.h"
#include "itkMagnitudeAndPhaseToComplexImageFilter.h"

#include "Util.h"

using namespace std;

enum class InputType { RealImag, MagPhase, Complex, Invalid };
template<typename TPixel> void Run(const InputType input,
                                   const string &r_name, const string &i_name,
                                   const string &m_name, const string &p_name,
                                   const string &x_name,
                                   const string &output, const vector<string> &o_names) {
	typedef itk::Image<TPixel, 4> TImage;
	typedef itk::Image<complex<TPixel>, 4> TXImage;

	typename TXImage::Pointer cImage = ITK_NULLPTR;
	switch (input) {
		case InputType::Complex: {
			auto read = itk::ImageFileReader<TXImage>::New();
			read->SetFileName(x_name);
			read->Update();
			cImage = read->GetOutput();
		} break;
		case InputType::RealImag: {
			auto read1 = itk::ImageFileReader<TImage>::New();
			read1->SetFileName(r_name);
			auto read2 = itk::ImageFileReader<TImage>::New();
			read2->SetFileName(i_name);
			auto compose = itk::ComposeImageFilter<TImage, TXImage>::New();
			compose->SetInput(0, read1->GetOutput());
			compose->SetInput(1, read2->GetOutput());
			compose->Update();
			cImage = compose->GetOutput();
		} break;
		case InputType::MagPhase: {
			auto read1 = itk::ImageFileReader<TImage>::New();
			auto read2 = itk::ImageFileReader<TImage>::New();
			read1->SetFileName(m_name);
			read2->SetFileName(p_name);
			auto compose = itk::MagnitudeAndPhaseToComplexImageFilter<TImage, TImage, TXImage>::New();
			compose->SetInput(0, read1->GetOutput());
			compose->SetInput(1, read2->GetOutput());
			compose->Update();
			cImage = compose->GetOutput();
		} break;
		case InputType::Invalid: throw(runtime_error("Should not happen.")); break;
	}

	auto scalarWriter = itk::ImageFileWriter<TImage>::New();
	for (size_t oi = 0; oi < output.size(); oi++) {
		char type = output.at(oi);
		scalarWriter->SetFileName(o_names.at(oi));
		switch (type) {
			case 'm': {
				auto o = itk::ComplexToModulusImageFilter<TXImage, TImage>::New();
				o->SetInput(cImage);
				scalarWriter->SetInput(o->GetOutput());
				scalarWriter->Update();
			} break;
			case 'p': {
				auto o = itk::ComplexToPhaseImageFilter<TXImage, TImage>::New();
				o->SetInput(cImage);
				scalarWriter->SetInput(o->GetOutput());
				scalarWriter->Update();
			} break;
			case 'r': {
				auto o = itk::ComplexToRealImageFilter<TXImage, TImage>::New();
				scalarWriter->SetInput(o->GetOutput());
				scalarWriter->Update();
			} break;
			case 'i': {
				auto o = itk::ComplexToImaginaryImageFilter<TXImage, TImage>::New();
				scalarWriter->SetInput(o->GetOutput());
				scalarWriter->Update();
			} break;
			case 'x': {
				auto complexWriter = itk::ImageFileWriter<TXImage>::New();
				complexWriter->SetFileName(o_names.at(oi));
				complexWriter->SetInput(cImage);
				complexWriter->Update();
			} break;
			default:
				throw(runtime_error("Unknown ouput type specifier: " + to_string(type)));
				break;
		}
	}
}


const string usage {
"Usage is: qcomplex [options] [input] -o [output] output_files \n\
\n\
Input must be specified as one of:\n\
	-r real_image -i imaginary_image\n\
	-m mag_image -p phase_image\n\
	-x complex_image\n\
Output can contain any/all of the letters 'mpric'. If you have more than one\n\
letter you must specify more than one output filename:\n\
	-o m : Output a magnitude image\n\
	   p : Output a phase image\n\
	   r : Output a real image\n\
	   i : Output an imaginary image\n\
	   x : Output a complex image\n\
Other options:\n\
	--dtype, -t f    : Datatype is float\n\
	           d     : Datatype is double\n\
	--fixge, -f      : Fix alternate slice, opposing phase issue on GE.\n\
Example:\n\
	qcomplex -m mag.nii -p phase.nii -o rix real.nii imag.nii complex.nii\n"
};

enum class PixelType { Float, Double };
const struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"dtype", required_argument, 0, 't'},
	{"fixge", no_argument, 0, 'f'},
	{0, 0, 0, 0}
};
const char* short_options = "hvr:i:m:p:x:o:t:f";

int main(int argc, char **argv) {
	bool verbose = false, forceDType = false, fixge = false;
	InputType inputType = InputType::Invalid;
	string output{""}, r_name, i_name, m_name, p_name, x_name;
	PixelType precision = PixelType::Float;
	InputType input = InputType::Invalid;
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
		case 'v': verbose = true; break;
		case 'r': r_name = string(optarg); break;
		case 'i': i_name = string(optarg); break;
		case 'm': m_name = string(optarg); break;
		case 'p': p_name = string(optarg); break;
		case 'x': x_name = string(optarg); break;
		case 'o': output = string(optarg); break;
		case 't':
			switch (*optarg) {
				case 'f': precision = PixelType::Float; break;
				case 'd': precision = PixelType::Double; break;
				default:
					cerr << "Unknown precision type " << optarg << endl;
					return EXIT_FAILURE;
			} break;
		case 'f': fixge = true; break;
		case 'h':
		case '?': // getopt will print an error message
			cout << usage << endl;
			return EXIT_SUCCESS;
		}
	}

	if ((argc - optind) != output.size()) {
		throw(runtime_error("Output type specifier does not match number of output filenames."));
	}
	vector<string> onames;
	while (optind < argc) {
		onames.push_back(string(argv[optind]));
		optind++;
	}

	if ((r_name.size() != 0) && (i_name.size() != 0)) {
		input = InputType::RealImag;
	} else if ((m_name.size() != 0) && (p_name.size() != 0)) {
		input = InputType::MagPhase;
	} else if (x_name.size() != 0) {
		input = InputType::Complex;
	} else {
		throw(runtime_error("Invalid combination of input specifiers."));
	}

	switch (precision) {
		case PixelType::Float: Run<float>(input, r_name, i_name, m_name, p_name, x_name, output, onames); break;
		case PixelType::Double: Run<double>(input, r_name, i_name, m_name, p_name, x_name, output, onames); break;
	}

	return EXIT_SUCCESS;
}

