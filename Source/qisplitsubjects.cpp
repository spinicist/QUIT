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
#include "itkShrinkImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkImageRegistrationMethod.h"

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

QI::ImageF::Pointer MaskWithLabel(const QI::ImageF::Pointer &image, const TLabelImage::Pointer &labels, const int l, const bool crop) {
    // convert the label image into a LabelMap
    typedef itk::LabelMap<itk::LabelObject<TLabel, 3>> TLabelMap;
    auto convert = itk::LabelImageToLabelMapFilter<TLabelImage, TLabelMap> ::New();
    convert->SetInput(labels);
    auto masker = itk::LabelMapMaskImageFilter<TLabelMap, QI::ImageF>::New();
    masker->SetInput(convert->GetOutput());
    masker->SetFeatureImage(image);
    masker->SetLabel(l);
    masker->SetBackgroundValue(0.);
    masker->SetNegated(false); // Mask outside the mask
    if (crop)
        masker->SetCrop(true);
    masker->Update();
    QI::ImageF::Pointer masked = masker->GetOutput();
    masked->DisconnectPipeline();
    return masked;
}

typedef itk::ImageMomentsCalculator<QI::ImageF> TMoments;
TMoments::VectorType GetCoG(const QI::ImageF::Pointer &img) {
    auto moments = TMoments::New();
    moments->SetImage(img);
    moments->Compute();
    return moments->GetCenterOfGravity();
}

typedef itk::Euler3DTransform<double> TRigid;
template<typename TImg, typename TInterpFunction>
typename TImg::Pointer ResampleImage(const typename TImg::Pointer &image, const typename TRigid::Pointer &tfm, const QI::ImageF::Pointer &reference = ITK_NULLPTR) {
    // Set these up to use in the loop
    typedef itk::ResampleImageFilter<TImg, TImg, double> TResampler;
    typedef TInterpFunction TInterpolator;
    typename TInterpolator::Pointer interp = TInterpolator::New();
    interp->SetInputImage(image);
    typename TResampler::Pointer resamp = TResampler::New();
    resamp->SetInput(image);
    resamp->SetInterpolator(interp);
    resamp->SetDefaultPixelValue(0.);
    resamp->SetTransform(tfm);
    if (reference)
        resamp->SetOutputParametersFromImage(reference);
    else
        resamp->SetOutputParametersFromImage(image);
    resamp->Update();
    typename TImg::Pointer rimage = resamp->GetOutput();
    rimage->DisconnectPipeline();
    return rimage;
}

void RegisterImageToReference(const QI::ImageF::Pointer &image, const QI::ImageF::Pointer &reference, TRigid::Pointer tfm, const int shrink) {
    typedef itk::SmoothingRecursiveGaussianImageFilter<QI::ImageF, QI::ImageF> TSmooth;
    typedef itk::ShrinkImageFilter<QI::ImageF, QI::ImageF> TShrink;
    typedef itk::RegularStepGradientDescentOptimizer TOpt;
    typedef itk::MattesMutualInformationImageToImageMetric<QI::ImageF, QI::ImageF> TMetric;
    typedef itk::ImageRegistrationMethod<QI::ImageF, QI::ImageF> TReg;
    typedef itk::LinearInterpolateImageFunction<QI::ImageF, double> TInterp;
    
    TSmooth::Pointer smooth_img = TSmooth::New();
    TSmooth::Pointer smooth_ref = TSmooth::New();
    TSmooth::SigmaArrayType sigma = image->GetSpacing();
    for (int i = 0; i < sigma.Size(); i++)
        sigma[i] = sigma[i] * shrink;
    smooth_img->SetInput(image);
    smooth_img->SetSigmaArray(sigma);
    smooth_ref->SetInput(reference);
    smooth_ref->SetSigmaArray(sigma);
    
    TShrink::Pointer shrink_img = TShrink::New();
    TShrink::Pointer shrink_ref = TShrink::New();
    shrink_img->SetInput(smooth_img->GetOutput());
    shrink_img->SetShrinkFactors(shrink);
    shrink_ref->SetInput(smooth_ref->GetOutput());
    shrink_ref->SetShrinkFactors(shrink);
    
    TMetric::Pointer metric = TMetric::New();
    TInterp::Pointer interp = TInterp::New();
    TOpt::Pointer opt = TOpt::New();
    TReg::Pointer reg = TReg::New();
    
    metric->SetNumberOfHistogramBins(32);
    metric->SetNumberOfSpatialSamples(10000);
    
    reg->SetMetric(metric);
    reg->SetOptimizer(opt);
    reg->SetTransform(tfm);
    reg->SetInterpolator(interp);
    reg->SetFixedImage(shrink_ref->GetOutput());
    reg->SetMovingImage(shrink_img->GetOutput());
    reg->SetFixedImageRegion(reference->GetLargestPossibleRegion());
    
    typedef TReg::ParametersType TPars;
    TPars initPars = tfm->GetParameters();
    
    TOpt::ScalesType scales(tfm->GetNumberOfParameters());
    const double translationScale = 1000.0;   // dynamic range of translations
    const double rotationScale    =    1.0;   // dynamic range of rotations
    scales[0] = 1.0 / rotationScale;
    scales[1] = 1.0 / rotationScale;
    scales[2] = 1.0 / rotationScale;
    scales[3] = 1.0 / translationScale;
    scales[4] = 1.0 / translationScale;
    scales[5] = 1.0 / translationScale;
    opt->SetScales(scales);
    opt->SetMaximumStepLength(4.0);
    opt->SetMinimumStepLength(0.01);
    opt->SetNumberOfIterations(250);
    
    reg->SetInitialTransformParameters(initPars);
    reg->Update();
    
    tfm->SetParameters(reg->GetLastTransformParameters());
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
Reference Options:\n\
    --ref, -r      : Specify a reference image for output space\n\
    --shrink, -S N : Specify shrink factor for registration step (default 2)\n\
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
    int keep = numeric_limits<int>::max(), size_threshold = 1000, shrink = 2, output_images = false;
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
        {"shrink", required_argument, 0, 'S'},
        {"ring", required_argument, 0, 'R'},
		{"rotX", required_argument, 0, 'X'},
		{"rotY", required_argument, 0, 'Y'},
		{"rotZ", required_argument, 0, 'Z'},
		{"oimgs", no_argument, &output_images, true},
		{0, 0, 0, 0}
	};
	const char* short_options = "hvr:c:k:s:S:";

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
        case 'k': keep = stoi(optarg); break;
        case 'S': shrink = stoi(optarg); break;
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

	for (auto i = 1; i <= keep; i++) {
        QI::ImageF::Pointer subject = MaskWithLabel(input, labels, i, true);
		TMoments::VectorType offset = -refCoG;
		double rotateAngle = 0.;
		if (alignment != ALIGN::NONE) {
            TMoments::VectorType CoG = GetCoG(subject);
			if (verbose) cout << "Subject " << i << " CoG is " << CoG << endl;
			rotateAngle = atan2(CoG[1], CoG[0]);
			if (alignment == ALIGN::RING_IN)
				rotateAngle = (M_PI / 2.) - rotateAngle;
			else if (alignment == ALIGN::RING_OUT)
				rotateAngle = (M_PI * 3./2.) - rotateAngle;
			if (verbose) cout << "Rotation angle is " << (rotateAngle*180./M_PI) << " degrees" << endl;
			offset += CoG;
		}

		auto tfm = TRigid::New();
		tfm->SetIdentity();
		tfm->SetRotation(angleX, angleY, angleZ - rotateAngle);
        tfm->SetOffset(offset);
        
        if (reference) {
            if (verbose) cout << "Registering to reference image..." << endl;
            RegisterImageToReference(subject, reference, tfm, shrink);
        }
        
        stringstream suffix; suffix << "_" << setfill('0') << setw(2) << i;
        fname = prefix + suffix.str() + ".tfm";
        if (verbose) cout << "Writing transform file " << fname << endl;
        auto tfmWriter = itk::TransformFileWriterTemplate<double>::New();
        tfmWriter->SetInput(tfm);
        tfmWriter->SetFileName(fname);
        tfmWriter->Update();
        
        if (output_images) {
            typedef itk::LinearInterpolateImageFunction<QI::ImageF, double> TLinInterp;
            typedef itk::NearestNeighborInterpolateImageFunction<TLabelImage, double> TNNInterp;
            QI::ImageF::Pointer rimage = ResampleImage<QI::ImageF, TLinInterp>(input, tfm, reference);
            TLabelImage::Pointer rlabels = ResampleImage<TLabelImage, TNNInterp>(labels, tfm, reference);
            QI::ImageF::Pointer masked = MaskWithLabel(rimage, rlabels, i, (reference == ITK_NULLPTR));
            fname = prefix + suffix.str() + ".nii";
            if (verbose) cout << "Writing output file " << fname << endl;
            QI::WriteImage(masked, fname);
		}
	}
	return EXIT_SUCCESS;
}

