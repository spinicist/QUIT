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
#include <sstream>
#include <iostream>
#include <iomanip>

#include "Util.h"
#include "Types.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageMomentsCalculator.h"
#include "itkAffineTransform.h"
#include "itkTransformFileWriter.h"

using namespace std;

int main(int argc, char **argv) {
	bool verbose = false;
	int indexptr = 0, c;
	int keep = 4, inwards = false, output_images = false;
	QI::ImageF::Pointer reference = ITK_NULLPTR;
	const string usage {
	"Usage is: qisplit input_file.nii [options]\n\
	\n\
	Options:\n\
		--help, -h    : Print this message\n\
		--verbose, -v : Print more information\n\
		--keep, -k    : Keep N largest objects (default 4)\n\
		--ref, -r     : Specify a reference image for CoG\n\
		--oimgs       : Output images\n\
		--inwards     : Subjects were scanned facing 'inwards'\n"
	};

	const struct option long_options[] =
	{
		{"help", no_argument, 0, 'h'},
		{"verbose", no_argument, 0, 'v'},
		{"keep", required_argument, 0, 'k'},
		{"ref", required_argument, 0, 'r'},
		{"oimgs", no_argument, &output_images, true},
		{"inwards", no_argument, &inwards, true},
		{0, 0, 0, 0}
	};
	const char* short_options = "hvr:k:";

	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
		case 'v': verbose = true; break;
		case 'h':
			cout << usage << endl;
			return EXIT_SUCCESS;
		case 'r': {
			auto refFile = QI::ReadImageF::New();
			refFile->SetFileName(optarg);
			refFile->Update();
			reference = refFile->GetOutput();
			reference->DisconnectPipeline();
		} break;
		case 'k':
			keep = atoi(optarg);
			break;
		case 0: // longopts flag
			break;
		case '?': // getopt will print an error message
			return EXIT_FAILURE;
		default:
			cout << "Unhandled option " << string(1, c) << endl;
			return EXIT_FAILURE;
		}
	}

	if ((argc - optind) != 1) {
		cout << usage << endl;
		throw(runtime_error("Wrong number of input arguments."));
	}

	auto input = QI::ReadImageF::New();
	string fname(argv[optind++]);
	input->SetFileName(fname);
	string prefix = QI::StripExt(fname);

	typedef unsigned int TLabel;
	typedef itk::Image<TLabel, 3> TLabelImage;
	/*typedef itk::BinaryThresholdImageFilter<QI::ImageF, TLabelImage> TThreshFilter;
	auto threshold = TThreshFilter::New();
	threshold->SetInput(input->GetOutput());
	threshold->SetLowerThreshold(atof(argv[optind++]));
	threshold->SetUpperThreshold(numeric_limits<float>::infinity());
	threshold->SetInsideValue(1);
	threshold->SetOutsideValue(0);*/
    // Use an Otsu Threshold filter to generate the mask
    if (verbose) cout << "Generating Otsu mask" << endl;
    auto otsuFilter = itk::OtsuThresholdImageFilter<QI::ImageF, TLabelImage>::New();
    otsuFilter->SetInput(input->GetOutput());
    otsuFilter->SetOutsideValue(1);
    otsuFilter->SetInsideValue(0);
    otsuFilter->Update();

	auto CC = itk::ConnectedComponentImageFilter<TLabelImage, TLabelImage>::New();
	CC->SetInput(otsuFilter->GetOutput());
	CC->Update();
	if (verbose) cout << "Found " << CC->GetObjectCount() << " objects in total, will keep " << keep << " largest." << endl;

	typedef itk::LabelShapeKeepNObjectsImageFilter<TLabelImage> TKeepN;
	auto keepN = TKeepN::New();
	keepN->SetInput(CC->GetOutput());
	keepN->SetBackgroundValue(0);
	keepN->SetNumberOfObjects(keep);
	keepN->SetAttribute(TKeepN::LabelObjectType::NUMBER_OF_PIXELS);

	auto relabel = itk::RelabelComponentImageFilter<TLabelImage, TLabelImage>::New();
	relabel->SetInput(keepN->GetOutput());

	if (verbose) cout << "Writing label image." << endl;
	auto output = itk::ImageFileWriter<TLabelImage>::New();
	output->SetInput(relabel->GetOutput());
	output->SetFileName(prefix + "_labels.nii");
	output->Update();

	auto labelStats = itk::LabelStatisticsImageFilter<QI::ImageF, TLabelImage>::New();
	labelStats->SetInput(input->GetOutput());
	labelStats->SetLabelInput(relabel->GetOutput());
	labelStats->Update();
	
	typedef itk::ImageMomentsCalculator<QI::ImageF> TMoments;
	TMoments::VectorType refCoG; refCoG.Fill(0);
	if (reference) {
		auto refMoments = TMoments::New();
		refMoments->SetImage(reference);
		refMoments->Compute();
		refCoG = refMoments->GetCenterOfGravity();
		if (verbose) cout << "Reference CoG is: " << refCoG << endl;
	}
	
	itk::Vector<double, 3> zAxis; zAxis.Fill(0); zAxis[2] = 1.0;
	for (auto i = 1; i <= 4; i++) {
		TLabelImage::RegionType region = labelStats->GetRegion(i);
		auto extract = itk::ExtractImageFilter<QI::ImageF, QI::ImageF>::New();
		extract->SetInput(input->GetOutput());
		extract->SetExtractionRegion(region);
		extract->SetDirectionCollapseToSubmatrix();
		extract->Update();

		auto moments = TMoments::New();
		moments->SetImage(extract->GetOutput());
		moments->Compute();
		TMoments::VectorType CoG = moments->GetCenterOfGravity();
		if (verbose) cout << "Writing object " << i << " CoG is " << CoG << endl;
		float rotateAngle = atan2(CoG[1], CoG[0]);
		if (inwards)
			rotateAngle = (M_PI / 2.) - rotateAngle;
		else
			rotateAngle = (M_PI * 3./2.) - rotateAngle;
		if (verbose) cout << "Rotation angle is " << (rotateAngle*180./M_PI) << " degrees" << endl;

		typedef itk::AffineTransform<double, 3> TAffine;
		auto tfm = TAffine::New();
		tfm->SetIdentity();
		tfm->SetOffset(CoG - refCoG);
		//tfm->SetCenter(CoG); // Want to leave this at 0 due to how ITK transforms work
		tfm->Rotate3D(zAxis, -rotateAngle, 1);

		if (verbose) cout << "Transform is: " << tfm << endl;

		stringstream suffix; suffix << "_" << setfill('0') << setw(2) << i;
		if (output_images) {
			fname = prefix + suffix.str() + ".nii";
			if (verbose) cout << "Writing output file " << fname << endl;
			auto output = itk::ImageFileWriter<QI::ImageF>::New();
			output->SetInput(extract->GetOutput());
			output->SetFileName(fname);
			output->Update();
		}

		fname = prefix + suffix.str() + ".tfm";
		if (verbose) cout << "Writing transform file " << fname << endl;
  		auto tfmWriter = itk::TransformFileWriterTemplate<double>::New();
  		tfmWriter->SetInput(tfm);
  		tfmWriter->SetFileName(fname);
  		tfmWriter->Update();
	}
	return EXIT_SUCCESS;
}

