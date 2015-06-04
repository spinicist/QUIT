/*
 *  qinewimg.cpp
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

#include "Util.h"

using namespace std;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: qinewimg filename dims voxdims [options]\n\
\n\
This is a tool to create Nifti files, either blank headers with orientation\n\
information, e.g. for registration, or files filled with simple patterns of\n\
data e.g. solid values, gradients, or blocks. The default is to create a\n\
3D file filled with zeros. Choose from the options below to create something\n\
else.\n\
\n\
Main Options:\n\
	--help, -h              : Print this message\n\
	--dims, -d \"N N N\"    : Set the dimensions (default 16 16 16)\n\
	--voxdims, -v \"X X X\" : Set the voxel dimensions (default 1mm iso)\n\
File Content Options:\n\
	--fill, -f X            : Fill the entire image with value X (default 0)\n\
	--grad, -g \"D L H\"    : Fill dimension D with a smooth gradient (low, high)\n\
	--step, -s \"D L H S\"  : Fill dimension D with stepped data (low, high, steps)\n"
};

/*
 * 	--precision, -p F       : Set the datatype to 32 bit float (default)\n\
	                I       : Set the datatype to 16 bit int\n\
	                C       : Set the datatype to 64 bit complex\n\
	--xform, -x FILE        : Copy header transform from another Nifti\n\
	--rank, -r N            : Set number of dimensions (max 5, default 3)\n\
*/

enum class FillTypes { Fill, Gradient, Steps };
enum class PixelTypes { Float, XFloat, Int16 };
//static bool verbose = false;
static int ndims = 3;
static FillTypes fillType = FillTypes::Fill;
static PixelTypes pixelType = PixelTypes::Float;
static float startVal = 0, deltaVal = 0, stopVal = 0;
static int stepLength = 1, fillDim = 0;
static struct option long_opts[] = {
	{"help",      no_argument,       0, 'h'},
	/*{"precision", required_argument, 0, 'p'},*/
	/*{"xform",     required_argument, 0, 'x'},*/
	/*{"rank",      required_argument, 0, 'r'},*/
	{"dims",      required_argument, 0, 'd'},
	{"voxdims",   required_argument, 0, 'v'},
	{"fill",      required_argument, 0, 'f'},
	{"grad",      required_argument, 0, 'g'},
	{"step",      required_argument, 0, 's'},
	{0, 0, 0, 0}
};
static const char* short_opts = "hd:v:f:g:s:"; /* p:x:r: */
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	QI::ImageF::Pointer newimg = QI::ImageF::New();
	QI::ImageF::RegionType imgRegion;
	QI::ImageF::IndexType imgIndex;
	QI::ImageF::SizeType imgSize;
	QI::ImageF::SpacingType imgSpacing;
	imgIndex.Fill(0);
	imgSize.Fill(0);
	imgSpacing.Fill(1.0);
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
			case 'p':
				/*switch (*optarg) {
					case 'F': pixelType = PixelTypes::Float; break;
					case 'C': pixelType = PixelTypes::XFloat; break;
					case 'I': pixelType = PixelTypes::Int16; break;
					default: throw(runtime_error("Unknown datatype: " + optarg)); break;
				}*/ break;
			case 'x': {
				//Nifti::File other(optarg);
				//xform = other.header().transform();
			} break;
			case 'r':
				/*ndims = atoi(optarg);
				if ((ndims < 2) || (ndims > 5)) {
					cerr << "Invalid number of dimensions. Must be 2-5." << endl;
					return EXIT_FAILURE;
				}*/
				break;
			case 'd': {
				string vals(optarg);
				istringstream stream(vals);
				for (int i = 0; i < 3; i++) {
					if (!(stream >> imgSize[i])) {
						throw(std::runtime_error("Failed to parse dims: " + vals));
					}
				}
			} break;
			case 'v': {
				string vals(optarg);
				istringstream stream(vals);
				for (int i = 0 ; i < 3; i++) {
					if (!(stream >> imgSpacing[i])) {
						throw(runtime_error("Failed to parse spacing: " + vals));
					}
				}
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
			case 's': {
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
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) < 1) {
		cerr << "Missing output filename." << endl;
		cout << usage << endl;
		return EXIT_FAILURE;
	} else if ((argc - optind) > 1) {
		cerr << "Unexpected extra arguments." << endl;
		cout << usage << endl;
		return EXIT_FAILURE;
	}
	string fName(argv[optind++]);

	imgRegion.SetIndex(imgIndex);
	imgRegion.SetSize(imgSize);
	newimg->SetRegions(imgRegion);
	newimg->SetSpacing(imgSpacing);
	newimg->Allocate();

	itk::ImageSliceIteratorWithIndex<QI::ImageF> it(newimg, imgRegion);

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
	QI::writeResult(newimg, fName);

	return EXIT_SUCCESS;
}


