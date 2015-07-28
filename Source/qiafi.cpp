/*
 *  qafi.cpp
 *
 *  Created by Tobias Wood on 2015/06/04.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <string>
#include <time.h>
#include <getopt.h>

#include "Types.h"
#include "Util.h"

#include "itkExtractImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"

using namespace std;

const string usage{
"Usage is: afi [options] input \n\
\
Options:\n\
	--verbose, -v   : Print more messages\n\
	--mask, -m file : Mask input with specified file.\n\
	--out, -o path  : Add a prefix to the output filenames.\n\
	--flip, -f      : Specify the nominal flip-angle (default 55 degrees)\n\
	--ratio, -r     : Specify TR2:TR1 ratio (default 5)\n"
};

static struct option long_options[] = {
	{"verbose", required_argument, 0, 'v'},
	{"mask",    required_argument, 0, 'm'},
	{"out",     required_argument, 0, 'o'},
	{"flip",    required_argument, 0, 'f'},
	{"ratio",   required_argument, 0, 'r'},
	{0, 0, 0, 0}
};
static const char *short_options = "vm:o:f:r:T:h";

template<class TPixel> class AFI {
public:
	AFI() {}
	~AFI() {}
	bool operator!=(const AFI &) const { return false; }
	bool operator==(const AFI &other) const { return !(*this != other); }

	inline TPixel operator()(const TPixel &r,
							 const TPixel &n) const
	{
		TPixel temp = (r*n - 1.) / (n - r);
		if (temp > 1.)
			temp = 1.;
		if (temp < -1.)
			temp = -1.;
		TPixel alpha = acos(temp) * 180. / M_PI;
		return alpha;
	}
};


int main(int argc, char **argv) {
	int indexptr = 0, c;
	string outPrefix = "AFI_";
	bool verbose = false;
	double n = 5., nomFlip = 55.;
	QI::ReadImageF::Pointer maskFile = ITK_NULLPTR;
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'm':
				cout << "Reading mask." << endl;
				maskFile = QI::ReadImageF::New();
				maskFile->SetFileName(optarg);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'f': nomFlip = atof(optarg); break;
			case 'r': n = atof(optarg); break;
			case 'T': itk::MultiThreader::SetGlobalDefaultNumberOfThreads(atoi(optarg)); break;
			case 'h':
				cout << usage << endl;
				return EXIT_SUCCESS;
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
			default:
				cout << "Unhandled option " << string(1, c) << endl;
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 1) {
		cout << usage << endl;
		return EXIT_FAILURE;
	}
	if (verbose) cout << "Opening input file " << argv[optind] << endl;
	auto inFile = QI::ReadTimeseriesF::New();
	inFile->SetFileName(argv[optind]);
	inFile->Update();
	if (verbose) {
		cout << "Nominal flip-angle is " << nomFlip << " degrees." << endl;
		cout << "TR2:TR1 ratio is " << n << endl;
	}
	auto volume1 = itk::ExtractImageFilter<QI::TimeseriesF, QI::ImageF>::New();
	auto volume2 = itk::ExtractImageFilter<QI::TimeseriesF, QI::ImageF>::New();
	auto region = inFile->GetOutput()->GetLargestPossibleRegion();
	region.GetModifiableSize()[3] = 0;
	region.GetModifiableIndex()[3] = 0;
	volume1->SetExtractionRegion(region);
	volume1->SetInput(inFile->GetOutput());
	volume1->SetDirectionCollapseToSubmatrix();
	region.GetModifiableIndex()[3] = 1;
	volume2->SetExtractionRegion(region);
	volume2->SetInput(inFile->GetOutput());
	volume2->SetDirectionCollapseToSubmatrix();

	auto imageRatio = itk::DivideImageFilter<QI::ImageF, QI::ImageF, QI::ImageF>::New();
	imageRatio->SetInput(0, volume2->GetOutput());
	imageRatio->SetInput(1, volume1->GetOutput());
	auto afi = itk::BinaryFunctorImageFilter<QI::ImageF, QI::ImageF, QI::ImageF, AFI<float>>::New();
	afi->SetInput1(imageRatio->GetOutput());
	afi->SetConstant2(n);
	auto B1 = itk::DivideImageFilter<QI::ImageF, QI::ImageF, QI::ImageF>::New();
	B1->SetInput1(afi->GetOutput());
	B1->SetConstant2(nomFlip);
	B1->Update();
    QI::writeResult(afi->GetOutput(), outPrefix + "angle" + QI::OutExt());
    QI::writeResult(B1->GetOutput(),  outPrefix + "B1" + QI::OutExt());
	if (verbose) cout << "Finished." << endl;
	return EXIT_SUCCESS;
}

