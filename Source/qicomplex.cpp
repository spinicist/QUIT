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
#include "Types.h"

using namespace std;

const string usage {
"Usage is: qcomplex [input options] [output options] [other options] \n\
\n\
Input is specified with lower case letters. One of the following\n\
combinations must be specified:\n\
	-m mag_image -p phase_image\n\
	-r real_image -i imaginary_image\n\
	-x complex_image\n\
\n\
Output is specified with upper case letters. One or more of the\n\
following can be specified:\n\
	-M : Output a magnitude image\n\
	-P : Output a phase image\n\
	-R : Output a real image\n\
	-I : Output an imaginary image\n\
	-X : Output a complex image\n\
\n\
Other options:\n\
	--double : Use double precision instead of float\n\
	--fixge  : Fix alternate slice problem with GE data.\n\
\n\
Example:\n\
	qicomplex -m mag.nii -p phase.nii -R real.nii -I imag.nii\n"
};
int verbose = false, use_double = false, fixge = false;
const struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"double", no_argument, &use_double, 1},
	{"fixge", no_argument, &fixge, 1},
	{0, 0, 0, 0}
};
const char* short_options = "hvm:M:p:P:r:R:i:I:x:X:";
int c, index_ptr = 0;

template<typename TPixel> void Run(int argc, char **argv) {
	typedef itk::Image<TPixel, 4>          TImage;
	typedef itk::Image<complex<TPixel>, 4> TXImage;
	typedef itk::ImageFileReader<TImage>   TReader;
	typedef itk::ImageFileReader<TXImage>  TXReader;
	typedef itk::ImageFileWriter<TImage>   TWriter;
	typedef itk::ImageFileWriter<TXImage>  TXWriter;

	typename TImage::Pointer img1 = ITK_NULLPTR, img2 = ITK_NULLPTR;
	typename TXImage::Pointer imgX = ITK_NULLPTR;

	if (verbose) cout << "Reading input files" << endl;
	bool ri = false, x = false;
	optind = 1;
	while ((c = getopt(argc, argv, short_options)) != -1) {
		switch (c) {
			case 'r':
				ri = true;
			case 'm': {
				auto read = TReader::New();
				read->SetFileName(optarg);
				read->Update();
				img1 = read->GetOutput();
			} break;
			case 'i':
				ri = true;
			case 'p': {
				auto read = TReader::New();
				read->SetFileName(optarg);
				read->Update();
				img2 = read->GetOutput();
			} break;
			case 'x': {
				x = true;
				typename TXReader::Pointer readX = TXReader::New();
				readX->SetFileName(optarg);
				readX->Update();
				imgX = readX->GetOutput();
			} break;
			default: break;
		}
	}

	if (x) {
		// Nothing to see here
	} else if (ri) {
		if (!(img1 && img2)) {
			throw(runtime_error("Must set real and imaginary inputs"));
		}
		if (verbose) cout << "Combining real and imaginary input" << endl;
		auto compose = itk::ComposeImageFilter<TImage, TXImage>::New();
		compose->SetInput(0, img1);
		compose->SetInput(1, img2);
		compose->Update();
		imgX = compose->GetOutput();
	} else {
		if (!(img1 && img2)) {
			throw(runtime_error("Must set magnitude and phase inputs"));
		}
		if (verbose) cout << "Combining magnitude and phase input" << endl;
		auto compose = itk::MagnitudeAndPhaseToComplexImageFilter<TImage, TImage, TXImage>::New();
		compose->SetInput(0, img1);
		compose->SetInput(1, img2);
		compose->Update();
		imgX = compose->GetOutput();
	}

	if (verbose) cout << "Writing output files" << endl;
	typename TWriter::Pointer write = TWriter::New();
	optind = 1;
	while ((c = getopt(argc, argv, short_options)) != -1) {
		switch (c) {
			case 'M': {
				auto o = itk::ComplexToModulusImageFilter<TXImage, TImage>::New();
				o->SetInput(imgX);
				write->SetFileName(optarg);
				write->SetInput(o->GetOutput());
				write->Update();
				if (verbose) cout << "Wrote magnitude image " + string(optarg) << endl;
			} break;
			case 'P': {
				auto o = itk::ComplexToPhaseImageFilter<TXImage, TImage>::New();
				o->SetInput(imgX);
				write->SetFileName(optarg);
				write->SetInput(o->GetOutput());
				write->Update();
				if (verbose) cout << "Wrote phase image " + string(optarg) << endl;
			} break;
			case 'R': {
				auto o = itk::ComplexToRealImageFilter<TXImage, TImage>::New();
				o->SetInput(imgX);
				write->SetFileName(optarg);
				write->SetInput(o->GetOutput());
				write->Update();
				if (verbose) cout << "Wrote real image " + string(optarg) << endl;
			} break;
			case 'I': {
				auto o = itk::ComplexToImaginaryImageFilter<TXImage, TImage>::New();
				o->SetInput(imgX);
				write->SetFileName(optarg);
				write->SetInput(o->GetOutput());
				write->Update();
				if (verbose) cout << "Wrote imaginary image " + string(optarg) << endl;
			} break;
			case 'X': {
				auto writeX = TXWriter::New();
				writeX->SetFileName(optarg);
				writeX->SetInput(imgX);
				writeX->Update();
				if (verbose) cout << "Wrote complex image " + string(optarg) << endl;
			} break;
			default: break;
		}
	}
}

int main(int argc, char **argv) {
	// Do one pass for the general options
	bool have_some_options = false;
	while ((c = getopt_long(argc, argv, short_options, long_options, &index_ptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'h':
				cout << usage << endl;
				return EXIT_SUCCESS;
			case 0: break; // A flag
			case 'm': case 'M': case 'p': case 'P':
			case 'r': case 'R': case 'i': case 'I':
			case 'x': case 'X':
				have_some_options = true;
				break;
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
			default:
			cout << "Unhandled option " << string(1, c) << endl;
			return EXIT_FAILURE;
		}
	}

	if (!have_some_options) {
		cout << usage << endl;
		return EXIT_FAILURE;
	}

	if (use_double) {
		if (verbose) cout << "Using double precision" << endl;
		Run<double>(argc, argv);
	} else {
		if (verbose) cout << "Using float precision" << endl;
		Run<float>(argc, argv);
	}

	return EXIT_SUCCESS;
}
