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
#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelMapMaskImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageMomentsCalculator.h"
#include "itkEuler3DTransform.h"
#include "itkTransformFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

using namespace std;

const string usage {
"Usage is: qisplit input_file.nii [options]\n\
\n\
Options:\n\
	--help, -h      : Print this message\n\
	--verbose, -v   : Print more information\n\
	--keep, -k N    : Keep N largest subjects\n\
    --size, -s N    : Only keep subjects with >N voxels (default 1000)\n\
	--ref, -r       : Specify a reference image for CoG\n\
	--center=IN/OUT : Center and rotate subjects that were scanned in standard\n\
	                  ring arrangement facing IN or OUT.\n\
	--rotX N        : Rotate by N degrees around the X axis.\n\
	--rotY N        : Rotate by N degrees around the Y axis.\n\
	--rotZ N        : Rotate by N degrees around the Z axis.\n\
	--oimgs         : Output images\n"
};
enum { CENTER_IN = 0, CENTER_OUT, CENTER_NONE};

int main(int argc, char **argv) {
	bool verbose = false;
	int indexptr = 0, c;
	int keep = numeric_limits<int>::max(), size_threshold = 1000, output_images = false,
	    center = CENTER_NONE;
    double angleX = 0., angleY = 0., angleZ = 0.;

	QI::ImageF::Pointer reference = ITK_NULLPTR;
	
	const struct option long_options[] =
	{
		{"help", no_argument, 0, 'h'},
		{"verbose", no_argument, 0, 'v'},
		{"keep", required_argument, 0, 'k'},
        {"size", required_argument, 0, 's'},
		{"ref", required_argument, 0, 'r'},
		{"center", required_argument, 0, 'c'},
		{"rotX", required_argument, 0, 'X'},
		{"rotY", required_argument, 0, 'Y'},
		{"rotZ", required_argument, 0, 'Z'},
		{"oimgs", no_argument, &output_images, true},
		{0, 0, 0, 0}
	};
	const char* short_options = "hvr:c:k:s:";

	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
		case 'v': verbose = true; break;
		case 'h':
			cout << usage << endl;
			return EXIT_SUCCESS;
		case 'r': reference = QI::ReadImage(optarg); break;
		case 'c':
			if (string(optarg) == "IN") {
				center = CENTER_IN;
			} else if (string(optarg) == "OUT") {
				center = CENTER_OUT;
			} else {
				cerr << "Unrecognised center specifier: " << string(optarg) << endl;
				return EXIT_FAILURE;
			} break;
		case 'k': keep = atoi(optarg); break;
        case 's': size_threshold = atoi(optarg); break;
        case 'X': angleX = stod(optarg)*M_PI/180.; break;
        case 'Y': angleY = stod(optarg)*M_PI/180.; break;
        case 'Z': angleZ = stod(optarg)*M_PI/180.; break;
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

	auto input = QI::ImageReaderF::New();
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
	
    typedef itk::RelabelComponentImageFilter<TLabelImage, TLabelImage> TRelabel;
    auto relabel = TRelabel::New();
    relabel->SetInput(CC->GetOutput());
    relabel->Update();
    // Relabel sorts on size by default, so now work out how many make the size threshold
    auto label_sizes = relabel->GetSizeOfObjectsInPixels();
    if (keep > label_sizes.size())
        keep = label_sizes.size();
    for (int i = 0; i < keep; i++) {
        if (label_sizes[i] < size_threshold) {
            keep = i;
            break;
        }
    }
    if (keep == 0) {
        cerr << "Found 0 objects to keep. Exiting." << endl;
        return EXIT_FAILURE;
    }
    if (verbose) cout << "Found " << keep << " objects to keep." << endl;

	typedef itk::LabelShapeKeepNObjectsImageFilter<TLabelImage> TKeepN;
	auto keepN = TKeepN::New();
	keepN->SetInput(relabel->GetOutput());
	keepN->SetBackgroundValue(0);
	keepN->SetNumberOfObjects(keep);
	keepN->SetAttribute(TKeepN::LabelObjectType::NUMBER_OF_PIXELS);

	if (verbose) cout << "Writing label image." << endl;
	auto output = itk::ImageFileWriter<TLabelImage>::New();
	output->SetInput(keepN->GetOutput());
	output->SetFileName(prefix + "_labels.nii");
	output->Update();

	typedef itk::LabelStatisticsImageFilter<QI::ImageF, TLabelImage> TLabelStats;
	auto labelStats = TLabelStats::New();
	labelStats->SetInput(input->GetOutput());
	labelStats->SetLabelInput(keepN->GetOutput());
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

    // Set these up to use in the loop
    typedef itk::ResampleImageFilter<TLabelImage, TLabelImage, double> TLabelResampler;
    typedef itk::NearestNeighborInterpolateImageFunction<TLabelImage, double> TLabelInterp;
    auto ilabels = TLabelInterp::New();
    ilabels->SetInputImage(keepN->GetOutput());
    auto rlabels = TLabelResampler::New();
    rlabels->SetInput(keepN->GetOutput());
    rlabels->SetInterpolator(ilabels);
    rlabels->SetDefaultPixelValue(0.);
    if (reference)
    	rlabels->SetOutputParametersFromImage(reference);
   	else
    	rlabels->SetOutputParametersFromImage(relabel->GetOutput());

    typedef itk::ResampleImageFilter<QI::ImageF, QI::ImageF, double> TResampleF;
    typedef itk::LinearInterpolateImageFunction<QI::ImageF, double> TLinearInterpF;
    auto interp = TLinearInterpF::New();
    interp->SetInputImage(input->GetOutput());
    auto rimage = TResampleF::New();
    rimage->SetInput(input->GetOutput());
    rimage->SetInterpolator(interp);
    rimage->SetDefaultPixelValue(0.);
    if (reference)
    	rimage->SetOutputParametersFromImage(reference);
    else
    	rimage->SetOutputParametersFromImage(input->GetOutput());

	for (auto i = 1; i <= keep; i++) {
		TLabelImage::RegionType region = labelStats->GetRegion(i);
		typedef itk::ExtractImageFilter<QI::ImageF, QI::ImageF> TExtractF;
		auto extract = TExtractF::New();
		extract->SetInput(input->GetOutput());
		extract->SetExtractionRegion(region);
		extract->SetDirectionCollapseToSubmatrix();
		extract->Update();

		TMoments::VectorType offset = -refCoG;
		double rotateAngle = 0.;
		if (center != CENTER_NONE) {
			auto moments = TMoments::New();
			moments->SetImage(extract->GetOutput());
			moments->Compute();
			TMoments::VectorType CoG = moments->GetCenterOfGravity();
			if (verbose) cout << "Subject " << i << " CoG is " << CoG << endl;
			rotateAngle = atan2(CoG[1], CoG[0]);
			if (center == CENTER_IN)
				rotateAngle = (M_PI / 2.) - rotateAngle;
			else if (center == CENTER_OUT)
				rotateAngle = (M_PI * 3./2.) - rotateAngle;
			if (verbose) cout << "Rotation angle is " << (rotateAngle*180./M_PI) << " degrees" << endl;
			offset += CoG;
		}

		typedef itk::Euler3DTransform<double> TRigid;
		auto tfm = TRigid::New();
		tfm->SetIdentity();
		tfm->SetOffset(offset);
		tfm->SetRotation(angleX, angleY, angleZ - rotateAngle);

		stringstream suffix; suffix << "_" << setfill('0') << setw(2) << i;
		if (output_images) {
            rlabels->SetTransform(tfm.GetPointer());
            rimage->SetTransform(tfm.GetPointer());
			
            rlabels->Update();
            rimage->Update();

              // convert the label image into a LabelMap
            typedef itk::LabelMap<itk::LabelObject<TLabel, 3>> TLabelMap;
            auto convert = itk::LabelImageToLabelMapFilter<TLabelImage, TLabelMap> ::New();
            convert->SetInput(rlabels->GetOutput());

            auto masker = itk::LabelMapMaskImageFilter<TLabelMap, QI::ImageF>::New();
            masker->SetInput(convert->GetOutput());
            masker->SetFeatureImage(rimage->GetOutput());
            masker->SetLabel(i);
            masker->SetBackgroundValue(0.);
            masker->SetNegated(false); // Mask outside the mask
            if (reference == ITK_NULLPTR)
            	masker->SetCrop(true);

            fname = prefix + suffix.str() + ".nii";
            if (verbose) cout << "Writing output file " << fname << endl;
            auto routput = itk::ImageFileWriter<QI::ImageF>::New();
            routput->SetFileName(fname);
            routput->SetInput(masker->GetOutput());
            routput->Update();
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

