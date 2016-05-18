/*
 *  qinewimage.cpp
 *
 *  Created by Tobias Wood on 02/06/2015.
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
#include "itkImageSliceIteratorWithIndex.h"

#include "QI/Util.h"

using namespace std;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: qinewimage filename dims voxdims [options]\n\
\n\
This is a tool to create 3D Nifti files, either blank headers with orientation\n\
information, e.g. for registration, or files filled with simple patterns of\n\
data e.g. solid values, gradients, or blocks. The default is to create a\n\
3D file filled with zeros. Choose from the options below to create something\n\
else.\n\
\n\
Main Options:\n\
	--help, -h              : Print this message\n\
	--size, -s \"N N N\"    : Set the dimensions (default 1 1 1)\n\
	--spacing, -p \"X Y Z\" : Set the voxel dimensions (default 1mm iso)\n\
	--origin, -o \"X Y Z\"  : Set the origin (default 0 0 0)\n\
File Content Options:\n\
	--fill, -f X            : Fill the entire image with value X (default 0)\n\
	--grad, -g \"D L H\"    : Fill dimension D with a smooth gradient (low, high)\n\
	--step, -s \"D L H S\"  : Fill dimension D with stepped data (low, high, steps)\n"
};

enum class FillTypes { Fill, Gradient, Steps };
int ndims = 3;
FillTypes fillType = FillTypes::Fill;
float startVal = 0, deltaVal = 0, stopVal = 0;
int stepLength = 1, fillDim = 0;
struct option long_opts[] = {
	{"help",      no_argument,       0, 'h'},
	{"size",      required_argument, 0, 's'},
	{"spacing",   required_argument, 0, 'p'},
	{"origin",    required_argument, 0, 'o'},
	{"fill",      required_argument, 0, 'f'},
	{"grad",      required_argument, 0, 'g'},
	{"step",      required_argument, 0, 't'},
	{0, 0, 0, 0}
};
static const char* short_opts = "hs:p:o:f:g:t:";
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	QI::VolumeF::Pointer newimg = QI::VolumeF::New();
	QI::VolumeF::RegionType  imgRegion;
	QI::VolumeF::IndexType   imgIndex;
	QI::VolumeF::SizeType    imgSize;
	QI::VolumeF::SpacingType imgSpacing;
	QI::VolumeF::PointType   imgOrigin;
	imgIndex.Fill(0);
	imgSize.Fill(1);
	imgSpacing.Fill(1.0);
	imgOrigin.Fill(0.);
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
			case 's': {
				string vals(optarg);
				istringstream stream(vals);
				for (int i = 0; i < 3; i++) {
					if (!(stream >> imgSize[i])) {
                        QI_EXCEPTION("Failed to parse dims: " << vals);
					}
				}
			} break;
			case 'p': {
				string vals(optarg);
				istringstream stream(vals);
				stream >> imgSpacing;
			} break;
			case 'o': {
				string vals(optarg);
				istringstream stream(vals);
				stream >> imgOrigin;
			} break;
			case 'f': fillType = FillTypes::Fill; startVal = atof(optarg); break;
			case 'g': {
				fillType = FillTypes::Gradient;
				string vals(optarg);
				istringstream stream(vals);
				stream >> fillDim;
				stream >> startVal;
				stream >> stopVal;
				deltaVal = (stopVal - startVal) / (imgSize[fillDim] - 1);
				stepLength = 1;
			} break;
			case 't': {
				fillType = FillTypes::Steps;
				string vals(optarg);
				istringstream stream(vals);
				stream >> fillDim;
				stream >> startVal;
				stream >> stopVal;
				int steps;
				stream >> steps;
				stepLength = imgSize[fillDim] / steps;
				deltaVal = (stopVal - startVal) / (steps - 1);
			} break;
			case 'h':
				cout << QI::GetVersion() << endl << usage << endl;
				return EXIT_SUCCESS;
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
			default:
				cout << "Unhandled option " << string(1, c) << endl;
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) < 1) {
		cerr << "Missing output filename." << endl;
		cout << QI::GetVersion() << endl << usage << endl;
		return EXIT_FAILURE;
	} else if ((argc - optind) > 1) {
		cerr << "Unexpected extra arguments." << endl;
		cout << QI::GetVersion() << endl << usage << endl;
		return EXIT_FAILURE;
	}
	string fName(argv[optind++]);

	imgRegion.SetIndex(imgIndex);
	imgRegion.SetSize(imgSize);
	newimg->SetRegions(imgRegion);
	newimg->SetSpacing(imgSpacing);
	newimg->SetOrigin(imgOrigin);
	newimg->Allocate();

	itk::ImageSliceIteratorWithIndex<QI::VolumeF> it(newimg, imgRegion);

	switch (fillDim) {
		case 0: it.SetFirstDirection(1); it.SetSecondDirection(2); break;
		case 1: it.SetFirstDirection(0); it.SetSecondDirection(2); break;
		case 2: it.SetFirstDirection(0); it.SetSecondDirection(1); break;
	}
	float val = startVal;
	it.GoToBegin();
	while (!it.IsAtEnd()) {
		while (!it.IsAtEndOfSlice()) {
			while (!it.IsAtEndOfLine()) {
				it.Set(val);
				++it;
			}
			it.NextLine();
		}
		it.NextSlice();
		if ((it.GetIndex()[fillDim] % stepLength) == (stepLength - 1)) val += deltaVal;
	}
    QI::WriteImage(newimg, fName);

	return EXIT_SUCCESS;
}


