/*
 *  qisplit.cpp
 *
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <getopt.h>
#include <string>
#include <iostream>

#include "Util.h"
#include "Types.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkOrientedBoundingBoxImageLabelMapFilter.h"
#include "itkAttributeImageLabelObject.h"
#include "itkOrientedBoundingBoxLabelObject.h"
#include "itkShapeLabelMapFilter.h"
#include "itkStatisticsLabelObject.h"
#include "itkLabelImageToLabelMapFilter.h"

using namespace std;

const string usage {
"Usage is: qisplit input_file.nii T N [options]\n\
\n\
N is the number of samples in the image\n\
T is the threshold\n\
Options:\n\
	--help, -h    : Print this message\n\
	--verbose, -v : Print more information\n"
};

const struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{0, 0, 0, 0}
};
const char* short_options = "hv";

int main(int argc, char **argv) {
	bool verbose = false;
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
		case 'v': verbose = true; break;
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

	if ((argc - optind) != 3) {
		cout << usage << endl;
		throw(runtime_error("Wrong number of input arguments."));
	}

	auto input = QI::ReadImageF::New();
	input->SetFileName(argv[optind++]);

	typedef unsigned int TLabel;
	typedef itk::Image<TLabel, 3> TLabelImage;
	typedef itk::BinaryThresholdImageFilter<QI::ImageF, TLabelImage> TThreshFilter;
	auto threshold = TThreshFilter::New();
	threshold->SetInput(input->GetOutput());
	threshold->SetLowerThreshold(atof(argv[optind++]));
	threshold->SetUpperThreshold(numeric_limits<float>::infinity());
	threshold->SetInsideValue(1);
	threshold->SetOutsideValue(0);

	auto CC = itk::ConnectedComponentImageFilter<TLabelImage, TLabelImage>::New();
	CC->SetInput(threshold->GetOutput());
	CC->Update();
	std::cout << "Number of objects: " << CC->GetObjectCount() << std::endl;

	typedef itk::LabelShapeKeepNObjectsImageFilter<TLabelImage> TKeepNFilter;
	auto keepN = TKeepNFilter::New();
	keepN->SetInput(CC->GetOutput());
	keepN->SetBackgroundValue(0);
	keepN->SetNumberOfObjects(atoi(argv[optind++]));
	keepN->SetAttribute(TKeepNFilter::LabelObjectType::NUMBER_OF_PIXELS);

	typedef itk::RelabelComponentImageFilter<TLabelImage, TLabelImage> TRelabelFilter;
	auto relabel = TRelabelFilter::New();
	relabel->SetInput(keepN->GetOutput());

	typedef itk::OrientedBoundingBoxLabelObject<TLabel, 3> TOBBLabelObject;
	typedef itk::AttributeImageLabelObject<TLabel, 3, QI::ImageF, TOBBLabelObject> TLabelObject;

	typedef itk::LabelMap<TLabelObject> TLabelMap;
	typedef itk::LabelImageToLabelMapFilter<TLabelImage, TLabelMap> TLabelToMapFilter;
	TLabelToMapFilter::Pointer toLabelMap = TLabelToMapFilter::New();

	toLabelMap->SetInput(relabel->GetOutput());
	toLabelMap->Update();

	typedef itk::OrientedBoundingBoxImageLabelMapFilter<TLabelMap> TOBBILabelMapFilter;

	TOBBILabelMapFilter::Pointer toOBBILabelMap = TOBBILabelMapFilter::New();
	toOBBILabelMap->SetInput(toLabelMap->GetOutput());
	toOBBILabelMap->SetAttributeImageSpacing(input->GetOutput()->GetSpacing());
	toOBBILabelMap->SetFeatureImage(input->GetOutput());
	toOBBILabelMap->Update();

	auto output = itk::ImageFileWriter<TLabelImage>::New();
	output->SetInput(relabel->GetOutput());
	output->SetFileName("all.nii");
	output->Update();

	unsigned int nLabels = toLabelMap->GetOutput()->GetNumberOfLabelObjects();
	cout << "There are " << nLabels << " labels." << endl;
	//toLabelMap->GetOutput()->PrintLabelObjects();
	//toOBBILabelMap->GetOutput()->PrintLabelObjects();
	for (auto i = 1; i <= nLabels; i++) {
		const TLabelObject *labelObject = toOBBILabelMap->GetOutput()->GetLabelObject(i);
		auto output = itk::ImageFileWriter<QI::ImageF>::New();
		output->SetInput(labelObject->GetAttributeImage());
		output->SetFileName("label" + to_string(i) + ".nii");
		output->Update();

		std::cout << "Label: " << i << std::endl;
		std::cout << "\tBBox: " << labelObject->GetBoundingBox()
		          << "\tCentroid: " << labelObject->GetCentroid()
		          << std::endl;
	}
	return EXIT_SUCCESS;
}

