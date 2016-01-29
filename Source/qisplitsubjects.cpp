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
#include "itkGrayscaleFillholeImageFilter.h"
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

typedef unsigned int TLabel;
typedef itk::Image<TLabel, 3> TLabelImage;

typename TLabelImage::Pointer ThresholdMask(const QI::ImageF::Pointer &img, const float thresh) {
    typedef itk::BinaryThresholdImageFilter<QI::ImageF, TLabelImage> TThreshFilter;
    auto threshold = TThreshFilter::New();
    threshold->SetInput(img);
    threshold->SetLowerThreshold(thresh);
    threshold->SetUpperThreshold(numeric_limits<float>::infinity());
    threshold->SetInsideValue(1);
    threshold->SetOutsideValue(0);
    threshold->Update();
    typename TLabelImage::Pointer mask = threshold->GetOutput();
    mask->DisconnectPipeline();
    return mask;
}

typename TLabelImage::Pointer OtsuMask(const QI::ImageF::Pointer &img) {
    auto otsuFilter = itk::OtsuThresholdImageFilter<QI::ImageF, TLabelImage>::New();
    otsuFilter->SetInput(img);
    otsuFilter->SetOutsideValue(1);
    otsuFilter->SetInsideValue(0);
    otsuFilter->Update();
    typename TLabelImage::Pointer mask = otsuFilter->GetOutput();
    mask->DisconnectPipeline();
    return mask;
}

typename TLabelImage::Pointer FindLabels(const TLabelImage::Pointer &mask, const int size_threshold, int &keep) {
    auto CC = itk::ConnectedComponentImageFilter<TLabelImage, TLabelImage>::New();
    auto relabel = itk::RelabelComponentImageFilter<TLabelImage, TLabelImage>::New();
    CC->SetInput(mask);
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
        throw(runtime_error("No labels found in mask"));
    }

    typedef itk::LabelShapeKeepNObjectsImageFilter<TLabelImage> TKeepN;
    TKeepN::Pointer keepN = TKeepN::New();
    keepN->SetInput(relabel->GetOutput());
    keepN->SetBackgroundValue(0);
    keepN->SetNumberOfObjects(keep);
    keepN->SetAttribute(TKeepN::LabelObjectType::NUMBER_OF_PIXELS);
    
    typedef itk::GrayscaleFillholeImageFilter<TLabelImage, TLabelImage> TFill;
    TFill::Pointer fill = TFill::New();
    fill->SetInput(keepN->GetOutput());
    fill->Update();
    
    typename TLabelImage::Pointer labels = fill->GetOutput();
    labels->DisconnectPipeline();
    return labels;
}

typedef itk::ImageMomentsCalculator<QI::ImageF> TMoments;
TMoments::VectorType GetCoG(const QI::ImageF::Pointer &img) {
    auto moments = TMoments::New();
    moments->SetImage(img);
    moments->Compute();
    return moments->GetCenterOfGravity();
}

const string usage {
"Usage is: qisplit input_file.nii [options]\n\
\n\
Options:\n\
    --help, -h     : Print this message\n\
    --verbose, -v  : Print more information\n\
    --keep, -k N   : Keep N largest subjects\n\
    --size, -s N   : Only keep subjects with >N voxels (default 1000)\n\
    --oimgs        : Output images (default only transforms)\n\
    --ref, -r      : Specify a reference image for output space\n\
Masking options (default is generate a mask with Otsu's method):\n\
    --mask, -m F   : Read the mask from file F\n\
    --thresh, -t N : Generate a mask by thresholding input at intensity N\n\
Alignment/arrangement correction (default is no correction):\n\
    --ring=IN/OUT  : Center and rotate subjects that were scanned in standard\n\
                     ring arrangement facing IN or OUT.\n\
    --rotX N       : Rotate by N degrees around the X axis.\n\
    --rotY N       : Rotate by N degrees around the Y axis.\n\
    --rotZ N       : Rotate by N degrees around the Z axis.\n\
"
};
enum class ALIGN { NONE = 0, RING_IN, RING_OUT };

int main(int argc, char **argv) {
	bool verbose = false;
	int indexptr = 0, c;
    int keep = numeric_limits<int>::max(), size_threshold = 1000, output_images = false;
    ALIGN alignment = ALIGN::NONE;
    float intensity_threshold = 0;
    double angleX = 0., angleY = 0., angleZ = 0.;

	QI::ImageF::Pointer reference = ITK_NULLPTR;
    TLabelImage::Pointer mask = ITK_NULLPTR;
	
	const struct option long_options[] =
	{
		{"help", no_argument, 0, 'h'},
		{"verbose", no_argument, 0, 'v'},
		{"keep", required_argument, 0, 'k'},
        {"size", required_argument, 0, 's'},
        {"thresh", required_argument, 0, 't'},
        {"mask", required_argument, 0, 'm'},
		{"ref", required_argument, 0, 'r'},
        {"ring", required_argument, 0, 'R'},
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
        case 'm': mask = QI::ReadImage<TLabelImage>(optarg); break;
        case 'R':
            if (string(optarg) == "IN") {
                alignment = ALIGN::RING_IN;
            } else if (string(optarg) == "OUT") {
                alignment = ALIGN::RING_OUT;
            } else {
                cerr << "Unrecognised ring alignment specifier: " << string(optarg) << endl;
                return EXIT_FAILURE;
            } break;
		case 'k': keep = atoi(optarg); break;
        case 's': size_threshold = atoi(optarg); break;
        case 't': intensity_threshold = stof(optarg); break;
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

    string fname(argv[optind++]);
    QI::ImageF::Pointer input = QI::ReadImage(fname);
    string prefix = QI::StripExt(fname);

    if (mask == ITK_NULLPTR) {
        if (intensity_threshold == 0)
            mask = OtsuMask(input);
        else
            mask = ThresholdMask(input, intensity_threshold);
    }

    TLabelImage::Pointer labels = FindLabels(mask, size_threshold, keep);
    if (verbose) cout << "Found " << keep << " subjects, saving labels." << endl;
    QI::WriteImage<TLabelImage>(labels, prefix + "_labels.nii");

	typedef itk::LabelStatisticsImageFilter<QI::ImageF, TLabelImage> TLabelStats;
	auto labelStats = TLabelStats::New();
	labelStats->SetInput(input);
	labelStats->SetLabelInput(labels);
	labelStats->Update();

    TMoments::VectorType refCoG;
    if (reference)
        refCoG = GetCoG(reference);
    else
        refCoG.Fill(0);
    
    // Set these up to use in the loop
    typedef itk::ResampleImageFilter<TLabelImage, TLabelImage, double> TLabelResampler;
    typedef itk::NearestNeighborInterpolateImageFunction<TLabelImage, double> TLabelInterp;
    auto ilabels = TLabelInterp::New();
    ilabels->SetInputImage(labels);
    auto rlabels = TLabelResampler::New();
    rlabels->SetInput(labels);
    rlabels->SetInterpolator(ilabels);
    rlabels->SetDefaultPixelValue(0.);
    if (reference)
    	rlabels->SetOutputParametersFromImage(reference);
   	else
    	rlabels->SetOutputParametersFromImage(labels);

    typedef itk::ResampleImageFilter<QI::ImageF, QI::ImageF, double> TResampleF;
    typedef itk::LinearInterpolateImageFunction<QI::ImageF, double> TLinearInterpF;
    auto interp = TLinearInterpF::New();
    interp->SetInputImage(input);
    auto rimage = TResampleF::New();
    rimage->SetInput(input);
    rimage->SetInterpolator(interp);
    rimage->SetDefaultPixelValue(0.);
    if (reference)
    	rimage->SetOutputParametersFromImage(reference);
    else
    	rimage->SetOutputParametersFromImage(input);

	for (auto i = 1; i <= keep; i++) {
		TLabelImage::RegionType region = labelStats->GetRegion(i);
		typedef itk::ExtractImageFilter<QI::ImageF, QI::ImageF> TExtractF;
		auto extract = TExtractF::New();
		extract->SetInput(input);
		extract->SetExtractionRegion(region);
		extract->SetDirectionCollapseToSubmatrix();
		extract->Update();

		TMoments::VectorType offset = -refCoG;
		double rotateAngle = 0.;
		if (alignment != ALIGN::NONE) {
            TMoments::VectorType CoG = GetCoG(extract->GetOutput());
			if (verbose) cout << "Subject " << i << " CoG is " << CoG << endl;
			rotateAngle = atan2(CoG[1], CoG[0]);
			if (alignment == ALIGN::RING_IN)
				rotateAngle = (M_PI / 2.) - rotateAngle;
			else if (alignment == ALIGN::RING_OUT)
				rotateAngle = (M_PI * 3./2.) - rotateAngle;
			if (verbose) cout << "Rotation angle is " << (rotateAngle*180./M_PI) << " degrees" << endl;
			offset += CoG;
		}

		typedef itk::Euler3DTransform<double> TRigid;
		auto tfm = TRigid::New();
		tfm->SetIdentity();
		tfm->SetRotation(angleX, angleY, angleZ - rotateAngle);
        tfm->SetOffset(offset);
        
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
            QI::WriteImage(masker->GetOutput(), fname);
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

