/*
 *  despot1.cpp
 *  Part of quitk
 *
 *  Created by Tobias Wood on 12/01/2015.
 *  Copyright (c) 2011-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <time.h>
#include <getopt.h>
#include <iostream>
#include <atomic>
#include <Eigen/Dense>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImage.h"
#include "itkImageToImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkComposeImageFilter.h"

#include "Filters/ImageToVectorFilter.h"
#include "Filters/DESPOT1Filter.h"
#include "Model.h"
#include "Sequence.h"
#include "Algorithm.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: despot1 [options] spgr_input \n\
\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--no-prompt, -n   : Suppress input prompts\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--mask, -m file   : Mask input with specified file\n\
	--B1, -b file     : B1 Map file (ratio)\n\
	--algo, -a l      : LLS algorithm (default)\n\
	           w      : WLLS algorithm\n\
	           n      : NLLS (Levenberg-Marquardt)\n\
	--its, -i N       : Max iterations for WLLS (default 4)\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

typedef itk::Image<float, 3> FloatImage;
typedef itk::VectorImage<float, 3> FloatVectorImage;
typedef itk::DESPOT1Filter<FloatVectorImage, FloatImage> DESPOT1;

static bool verbose = false, prompt = true, all_residuals = false;
static size_t nIterations = 4;
static string outPrefix;
static shared_ptr<D1> algo = make_shared<D1>();
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"B1", required_argument, 0, 'b'},
	{"algo", required_argument, 0, 'a'},
	{"its", required_argument, 0, 'i'},
	{"threads", required_argument, 0, 'T'},
	{"resids", no_argument, 0, 'r'},
	{0, 0, 0, 0}
};
static const char *short_opts = "hvnm:o:b:a:i:T:r";
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	try { // To fix uncaught exceptions on Mac
	//cout << version << endl << credit_shared << endl;
	Eigen::initParallel();

	typedef itk::ImageFileReader<FloatImage> Reader;
	typedef itk::ImageFileReader<itk::Image<float, 4>> Reader4D;
	typedef itk::ImageFileWriter<FloatImage> Writer;

	Reader::Pointer mask = ITK_NULLPTR;
	Reader::Pointer B1   = ITK_NULLPTR;

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'm':
				cout << "Opening mask file " << optarg << endl;
				mask = Reader::New();
				mask->SetFileName(optarg);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'b':
				cout << "Opening B1 file: " << optarg << endl;
				B1 = Reader::New();
				B1->SetFileName(optarg);
				break;
			case 'a':
				switch (*optarg) {
					case 'l': algo->setType(D1::Type::LLS);  cout << "LLS algorithm selected." << endl; break;
					case 'w': algo->setType(D1::Type::WLLS); cout << "WLLS algorithm selected." << endl; break;
					case 'n': algo->setType(D1::Type::NLLS); cout << "NLLS algorithm selected." << endl; break;
					default:
						cout << "Unknown algorithm type " << optarg << endl;
						return EXIT_FAILURE;
						break;
				} break;
			case 'i':
				nIterations = atoi(optarg);
				break;
			case 'T':
				//threads.resize(atoi(optarg));
				break;
			case 'r': all_residuals = true; break;
			case 'h':
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 1) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		return EXIT_FAILURE;
	}

	//**************************************************************************
	#pragma mark Gather SPGR data
	//**************************************************************************
	cout << "Opening SPGR file: " << argv[optind] << endl;
	Reader4D::Pointer input = Reader4D::New();
	input->SetFileName(argv[optind]);

	typedef ImageToVectorFilter<float> Converter;
	cout << "Creating converter" << endl;
	Converter::Pointer convert = Converter::New();
	convert->SetInput(input->GetOutput());

	shared_ptr<SPGRSimple> spgrSequence = make_shared<SPGRSimple>(prompt);
	if (verbose) {
		cout << spgrSequence;
		cout << "Ouput prefix will be: " << outPrefix << endl;
	}

	auto vectorImage = convert->GetOutput();
	cout << "Creating DESPOT1 Filter" << endl;
	DESPOT1::Pointer d1 = DESPOT1::New();
	d1->SetInput(convert->GetOutput());
	if (mask)
		d1->SetMask(mask->GetOutput());
	if (B1)
		d1->SetB1(B1->GetOutput());
	d1->SetSequence(spgrSequence);
	algo->setIterations(nIterations);
	d1->SetAlgorithm(algo);
	cout << "Created filter" << endl;

	if (verbose)
		cout << "Writing results." << endl;
	outPrefix = outPrefix + "D1_";

	cout << "Creating Writers" << endl;
	Writer::Pointer T1File = Writer::New();
	Writer::Pointer PDFile = Writer::New();
	Writer::Pointer ResFile = Writer::New();

	T1File->SetFileName(outPrefix + "T1.nii");
	PDFile->SetFileName(outPrefix + "PD.nii");
	ResFile->SetFileName(outPrefix + "residual.nii");

	T1File->SetInput(d1->GetOutput(0));
	PDFile->SetInput(d1->GetOutput(1));
	ResFile->SetInput(d1->GetOutput(2));

	cout << "Processing" << endl;
	d1->Update();
	cout << "Writing output files" << endl;
	T1File->Update();
	PDFile->Update();
	ResFile->Update();

	/*
	if (all_residuals) {
		outHdr.intent_name = "Residuals";
		outHdr.setDim(4, spgrSequence.size());
		outFile.setHeader(outHdr);
		outFile.open(outPrefix + "residuals" + OutExt(), Nifti::Mode::Write);
		outFile.writeVolumes(ResidsVols.begin(), ResidsVols.end(), 0, spgrSequence.size());
		outFile.close();
	}*/
	cout << "All done." << endl;
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
