/*
 *  qisplitsubjects.cpp
 *
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <array>
#include <getopt.h>
#include <iomanip>
#include <iostream>

#include "ImageIO.h"
#include "ImageTypes.h"
#include "Masking.h"
#include "Util.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkConstantBoundaryCondition.h"
#include "itkEuler3DTransform.h"
#include "itkExtractImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkImageMomentsCalculator.h"
#include "itkImageRegistrationMethod.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelMapMaskImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkTransformFileWriter.h"
#include "itkWindowedSincInterpolateImageFunction.h"

using namespace std;

QI::VolumeF::Pointer MaskWithLabel(const QI::VolumeF::Pointer &image,
                                   const QI::VolumeI::Pointer &labels, const int l,
                                   const bool crop) {
    // convert the label image into a LabelMap
    typedef itk::LabelMap<itk::LabelObject<int, 3>> TLabelMap;
    auto convert = itk::LabelImageToLabelMapFilter<QI::VolumeI, TLabelMap>::New();
    convert->SetInput(labels);
    auto masker = itk::LabelMapMaskImageFilter<TLabelMap, QI::VolumeF>::New();
    masker->SetInput(convert->GetOutput());
    masker->SetFeatureImage(image);
    masker->SetLabel(l);
    masker->SetBackgroundValue(0.);
    masker->SetNegated(false); // Mask outside the mask
    if (crop)
        masker->SetCrop(true);
    masker->Update();
    QI::VolumeF::Pointer masked = masker->GetOutput();
    masked->DisconnectPipeline();
    return masked;
}

typedef itk::ImageMomentsCalculator<QI::VolumeF> TMoments;
TMoments::VectorType                             GetCoG(const QI::VolumeF::Pointer &img) {
    auto moments = TMoments::New();
    moments->SetImage(img);
    moments->Compute();
    return -moments->GetCenterOfGravity(); // ITK seems to put a negative sign on the CoG
}

typedef itk::Euler3DTransform<double> TRigid;
template <typename TImg, typename TInterp>
typename TImg::Pointer ResampleImage(const typename TImg::Pointer &  image,
                                     const typename TRigid::Pointer &tfm,
                                     const QI::VolumeF::Pointer &    reference = ITK_NULLPTR) {
    typedef itk::ResampleImageFilter<TImg, TImg, double> TResampler;
    typename TInterp::Pointer                            interp = TInterp::New();
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
    // Get rid of any negative values
    typedef itk::ThresholdImageFilter<TImg> TThreshold;
    auto                                    threshold = TThreshold::New();
    threshold->SetInput(resamp->GetOutput());
    threshold->ThresholdBelow(0);
    threshold->SetOutsideValue(0);
    threshold->Update();
    typename TImg::Pointer rimage = threshold->GetOutput();
    rimage->DisconnectPipeline();
    return rimage;
}

/*
 * Typedefs for registration
 */
typedef itk::SmoothingRecursiveGaussianImageFilter<QI::VolumeF, QI::VolumeF>     TSmooth;
typedef itk::ShrinkImageFilter<QI::VolumeF, QI::VolumeF>                         TShrink;
typedef itk::RegularStepGradientDescentOptimizer                                 TOpt;
typedef itk::MattesMutualInformationImageToImageMetric<QI::VolumeF, QI::VolumeF> TMetric;
typedef itk::ImageRegistrationMethod<QI::VolumeF, QI::VolumeF>                   TReg;
typedef TReg::ParametersType                                                     TPars;
typedef itk::LinearInterpolateImageFunction<QI::VolumeF, double>                 TInterp;

/*
 * Helper functions for registration
 */
TOpt::ScalesType MakeScales(double rotScale, double tScale) {
    TOpt::ScalesType scales(TRigid::ParametersDimension);
    scales[0] = rotScale;
    scales[1] = rotScale;
    scales[2] = rotScale;
    scales[3] = tScale;
    scales[4] = tScale;
    scales[5] = tScale;
    return scales;
}

TShrink::ShrinkFactorsType MakeShrink(const double &              gridSpacing,
                                      const QI::VolumeF::Pointer &image) {
    TShrink::ShrinkFactorsType shrink;
    for (size_t i = 0; i < shrink.Size(); i++) {
        shrink[i] = round(gridSpacing / image->GetSpacing()[i]);
        if (shrink[i] < 1)
            shrink[i] = 1;
    }
    return shrink;
}

TPars MakePars(const TPars &ip, double ax, double ay, double az, double tx, double ty, double tz) {
    TPars p((TRigid::ParametersDimension));
    p[0] = ip[0] + ax;
    p[1] = ip[1] + ay;
    p[2] = ip[2] + az;
    p[3] = ip[3] + tx;
    p[4] = ip[4] + ty;
    p[5] = ip[5] + tz;
    return p;
}

/*
 * Actual registration step
 */
void RegisterImageToReference(const QI::VolumeF::Pointer &image,
                              const QI::VolumeF::Pointer &reference, TRigid::Pointer tfm,
                              double gridSpacing, double angleStep, int searchAngles,
                              const int iterations, const bool verbose) {

    if (verbose)
        cout << "Rescaling to matched intensity ranges" << endl;
    typedef itk::RescaleIntensityImageFilter<QI::VolumeF, QI::VolumeF> TRescale;
    TRescale::Pointer scale_image = TRescale::New();
    TRescale::Pointer scale_ref   = TRescale::New();
    scale_image->SetInput(image);
    scale_ref->SetInput(reference);
    const double desiredMinimum = 0.0;
    const double desiredMaximum = 1000.0;
    scale_image->SetOutputMinimum(desiredMinimum);
    scale_image->SetOutputMaximum(desiredMaximum);
    scale_image->UpdateLargestPossibleRegion();
    scale_ref->SetOutputMinimum(desiredMinimum);
    scale_ref->SetOutputMaximum(desiredMaximum);
    scale_ref->UpdateLargestPossibleRegion();

    if (verbose)
        cout << "Matching histograms" << endl;
    typedef itk::HistogramMatchingImageFilter<QI::VolumeF, QI::VolumeF> THistMatch;
    THistMatch::Pointer hist_match = THistMatch::New();
    hist_match->SetReferenceImage(scale_ref->GetOutput());
    hist_match->SetInput(scale_image->GetOutput());
    hist_match->SetNumberOfHistogramLevels(100);
    hist_match->SetNumberOfMatchPoints(15);
    hist_match->ThresholdAtMeanIntensityOn();
    hist_match->Update();

    TSmooth::Pointer smooth_img = TSmooth::New();
    smooth_img->SetInput(hist_match->GetOutput());
    TSmooth::Pointer smooth_ref = TSmooth::New();
    smooth_ref->SetInput(scale_ref->GetOutput());
    TShrink::Pointer shrink_img = TShrink::New();
    shrink_img->SetInput(smooth_img->GetOutput());
    TShrink::Pointer shrink_ref = TShrink::New();
    shrink_ref->SetInput(smooth_ref->GetOutput());

    TMetric::Pointer metric = TMetric::New();
    metric->SetNumberOfHistogramBins(32);
    metric->SetNumberOfSpatialSamples(10000);
    TInterp::Pointer interp = TInterp::New();

    TOpt::Pointer opt = TOpt::New();
    opt->SetScales(MakeScales(1.0, 1. / 1000.));
    opt->SetMaximumStepLength(1.0);
    opt->SetMinimumStepLength(0.01);
    opt->SetNumberOfIterations(iterations);

    TReg::Pointer reg = TReg::New();
    reg->SetMetric(metric);
    reg->SetOptimizer(opt);
    reg->SetTransform(tfm);
    reg->SetInterpolator(interp);
    reg->SetFixedImage(shrink_ref->GetOutput());
    reg->SetMovingImage(shrink_img->GetOutput());
    reg->SetFixedImageRegion(reference->GetLargestPossibleRegion());

    double pangle     = angleStep * M_PI / 180.;
    double initMetric = 0, bestMetric = 0;
    TPars  initPars;
    TPars  bestPars = tfm->GetParameters();

    if (verbose)
        cout << "Starting registration" << endl;
    std::array<std::array<int, 3>, 7> translations{{{{0, 0, 0}},
                                                    {{-1, 0, 0}},
                                                    {{1, 0, 0}},
                                                    {{0, -1, 0}},
                                                    {{0, 1, 0}},
                                                    {{0, 0, -1}},
                                                    {{0, 0, 1}}}};
    do {
        TShrink::ShrinkFactorsType imageShrink = MakeShrink(gridSpacing, image);
        shrink_img->SetShrinkFactors(imageShrink);
        TShrink::ShrinkFactorsType refShrink = MakeShrink(gridSpacing, reference);
        shrink_ref->SetShrinkFactors(refShrink);

        TSmooth::SigmaArrayType smooth;
        smooth.Fill(gridSpacing / 2);
        smooth_img->SetSigmaArray(smooth);
        smooth_ref->SetSigmaArray(smooth);

        if (verbose)
            cout << "Grid: " << gridSpacing << " Image Shrink: " << imageShrink
                 << " Ref Shrink:   " << refShrink << endl;

        initPars = bestPars;
        reg->SetInitialTransformParameters(initPars);
        try {
            reg->Update();
            initMetric = opt->GetValue();
            bestMetric = initMetric;
            bestPars   = reg->GetLastTransformParameters();
            if (verbose)
                cout << "Initial metric at this level: " << initMetric << endl;
        } catch (itk::ExceptionObject &e) {
            if (verbose)
                cout << "Initial registration failed with parameters: "
                     << reg->GetLastTransformParameters() << endl;
        }
        for (double ax = -searchAngles * pangle; ax <= searchAngles * pangle; ax += pangle) {
            for (double ay = -searchAngles * pangle; ay <= searchAngles * pangle; ay += pangle) {
                for (double az = -searchAngles * pangle; az <= searchAngles * pangle;
                     az += pangle) {
                    for (const auto &t : translations) {
                        TPars p = MakePars(initPars, ax, ay, az, t[0] * gridSpacing,
                                           t[1] * gridSpacing, t[2] * gridSpacing);
                        reg->SetInitialTransformParameters(p);
                        try {
                            reg->Update();
                        } catch (itk::ExceptionObject &e) {
                            if (verbose)
                                cout << "Registration failed for parameters: "
                                     << reg->GetLastTransformParameters() << endl;
                        }
                        if (opt->GetValue() < bestMetric) {
                            bestMetric = opt->GetValue();
                            bestPars   = reg->GetLastTransformParameters();
                            if (verbose)
                                cout << "Metric improved to: " << bestMetric
                                     << " Iterations: " << opt->GetCurrentIteration() << endl;
                        }
                    }
                }
            }
        }
        pangle /= 2;
        gridSpacing /= 2;
    } while ((gridSpacing >= image->GetSpacing()[0]) && (bestMetric != initMetric));
    if (verbose)
        cout << "Finished" << endl;
    tfm->SetParameters(bestPars);
}

const string usage{"Usage is: qisplit input_file [options]\n\
\n\
Options:\n\
    --help, -h     : Print this message\n\
    --verbose, -v  : Print more information\n\
    --keep, -k N   : Keep N largest subjects\n\
    --size, -s N   : Only keep subjects with >N voxels (default 1000)\n\
    --oimgs        : Output images (default only transforms)\n\
Reference Options:\n\
    --ref, -R FILE : Specify a reference image for output space\n\
    --grid, -G N   : Specify initial grid scale (default 1mm)\n\
    --iters, -I N  : Specify the max number of iterations (default 25)\n\
    --angle, -A N  : Search angle step in degrees (default 30)\n\
    --nangle, -N N : Number of angle steps (default 1)\n\
Masking options (default is generate a mask with Otsu's method):\n\
    --mask, -m F   : Read the mask from file F\n\
    --thresh, -t N : Generate a mask by thresholding input at intensity N\n\
Alignment/arrangement correction (default is no correction):\n\
    --align=IN     : Ring arrangement facing IN\n\
            OUT    : Ring arrangement facing OUT\n\
            READ   : Read rotation angle from stdin\n\
    --rotX N       : Rotate by N degrees around the X axis.\n\
    --rotY N       : Rotate by N degrees around the Y axis.\n\
    --rotZ N       : Rotate by N degrees around the Z axis.\n\
"};
enum class ALIGN { NONE = 0, RING_IN, RING_OUT, READ };

int main(int argc, char **argv) {
    bool verbose  = false;
    int  indexptr = 0, c;
    int  keep = numeric_limits<int>::max(), size_threshold = 1000, output_images = false,
        iterations = 25, angleSteps = 1;
    ALIGN  alignment           = ALIGN::NONE;
    float  intensity_threshold = 0;
    double angleX = 0., angleY = 0., angleZ = 0., angleStep = 30., gridSpacing = 1.0;

    QI::VolumeF::Pointer reference = ITK_NULLPTR;
    QI::VolumeI::Pointer mask      = ITK_NULLPTR;

    const struct option long_options[] = {{"help", no_argument, 0, 'h'},
                                          {"verbose", no_argument, 0, 'v'},
                                          {"keep", required_argument, 0, 'k'},
                                          {"size", required_argument, 0, 's'},
                                          {"thresh", required_argument, 0, 't'},
                                          {"mask", required_argument, 0, 'm'},
                                          {"ref", required_argument, 0, 'R'},
                                          {"grid", required_argument, 0, 'G'},
                                          {"iters", required_argument, 0, 'I'},
                                          {"align", required_argument, 0, 'a'},
                                          {"angle", required_argument, 0, 'A'},
                                          {"nangle", required_argument, 0, 'N'},
                                          {"rotX", required_argument, 0, 'X'},
                                          {"rotY", required_argument, 0, 'Y'},
                                          {"rotZ", required_argument, 0, 'Z'},
                                          {"oimgs", no_argument, &output_images, true},
                                          {0, 0, 0, 0}};
    const char *        short_options  = "hvR:G:I:A:N:c:k:s:";

    while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
        switch (c) {
        case 'v':
            verbose = true;
            break;
        case 'h':
            cout << QI::GetVersion() << endl << usage << endl;
            return EXIT_SUCCESS;
        case 'm':
            mask = QI::ReadImage<QI::VolumeI>(optarg);
            break;
        case 'a':
            if (string(optarg) == "IN") {
                alignment = ALIGN::RING_IN;
            } else if (string(optarg) == "OUT") {
                alignment = ALIGN::RING_OUT;
            } else if (string(optarg) == "READ") {
                alignment = ALIGN::READ;
            } else {
                cerr << "Unrecognised alignment specifier: " << string(optarg) << endl;
                return EXIT_FAILURE;
            }
            break;
        case 'k':
            keep = stoi(optarg);
            break;
        case 'R':
            reference = QI::ReadImage(optarg);
            break;
        case 'G':
            gridSpacing = stod(optarg);
            break;
        case 'I':
            iterations = stoi(optarg);
            break;
        case 'A':
            angleStep = stod(optarg);
            break;
        case 'N':
            angleSteps = stoi(optarg);
            break;
        case 's':
            size_threshold = atoi(optarg);
            break;
        case 't':
            intensity_threshold = stof(optarg);
            break;
        case 'X':
            angleX = stod(optarg) * M_PI / 180.;
            break;
        case 'Y':
            angleY = stod(optarg) * M_PI / 180.;
            break;
        case 'Z':
            angleZ = stod(optarg) * M_PI / 180.;
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
        cout << QI::GetVersion() << endl << usage << endl;
        QI::Fail("Wrong number of input arguments.");
    }

    string               fname(argv[optind++]);
    QI::VolumeF::Pointer input  = QI::ReadImage(fname);
    string               prefix = QI::StripExt(fname);

    if (!mask) {
        if (intensity_threshold == 0)
            mask = QI::OtsuMask(input);
        else
            mask = QI::ThresholdMask(input, intensity_threshold);
    }

    QI::VolumeI::Pointer labels;
    QI::FindLabels(mask, size_threshold, keep, labels);
    if (verbose)
        cout << "Found " << keep << " subjects, saving labels." << endl;
    QI::WriteImage(labels, prefix + "_labels" + QI::OutExt());

    typedef itk::LabelStatisticsImageFilter<QI::VolumeF, QI::VolumeI> TLabelStats;
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
        QI::VolumeF::Pointer subject = MaskWithLabel(input, labels, i, true);
        TMoments::VectorType CoG     = GetCoG(subject);
        TMoments::VectorType offset  = refCoG - CoG;
        if (verbose)
            cout << "Subject " << i << " CoG is " << CoG << ", angle is "
                 << (180. / M_PI) * atan2(CoG[1], CoG[0]) << endl;
        double rotZ = 0.;
        switch (alignment) {
        case ALIGN::NONE:
            break;
        case ALIGN::RING_IN:
            rotZ = (M_PI / 2. + atan2(CoG[1], CoG[0]));
            break;
        case ALIGN::RING_OUT:
            rotZ = (-M_PI / 2. + atan2(CoG[1], CoG[0]));
            break;
        case ALIGN::READ:
            cin >> angleX;
            angleX *= (M_PI / 180.);
            cin >> angleY;
            angleY *= (M_PI / 180.);
            cin >> angleZ;
            angleZ *= (M_PI / 180.);
            break;
        }
        if (verbose)
            cout << "Initial rotation angles are " << (angleX * 180. / M_PI) << " "
                 << (angleY * 180. / M_PI) << " " << (angleZ * 180. / M_PI) << " degrees" << endl;

        TRigid::Pointer tfm = TRigid::New();
        tfm->SetIdentity();
        tfm->SetRotation(angleX, angleY, angleZ + rotZ);
        tfm->SetOffset(offset);

        if (reference) {
            if (verbose)
                cout << "Registering to reference image..." << endl;
            RegisterImageToReference(subject, reference, tfm, gridSpacing, angleStep, angleSteps,
                                     iterations, verbose);
        }

        stringstream suffix;
        suffix << "_" << setfill('0') << setw(2) << i;
        fname = prefix + suffix.str() + ".tfm";
        if (verbose)
            cout << "Writing transform file " << fname << endl;
        auto tfmWriter = itk::TransformFileWriterTemplate<double>::New();
        tfmWriter->SetInput(tfm);
        tfmWriter->SetFileName(fname);
        tfmWriter->Update();

        if (output_images) {
            typedef itk::WindowedSincInterpolateImageFunction<
                QI::VolumeF, 5, itk::Function::LanczosWindowFunction<5>,
                itk::ConstantBoundaryCondition<QI::VolumeF>, double>
                                                                                      TInterp;
            typedef itk::NearestNeighborInterpolateImageFunction<QI::VolumeI, double> TNNInterp;
            if (verbose)
                cout << "Resampling image" << endl;
            QI::VolumeF::Pointer rimage =
                ResampleImage<QI::VolumeF, TInterp>(subject, tfm, reference);
            QI::VolumeI::Pointer rlabels =
                ResampleImage<QI::VolumeI, TNNInterp>(labels, tfm, reference);
            typedef itk::BinaryThresholdImageFilter<QI::VolumeI, QI::VolumeI> TThreshFilter;
            auto rthresh = TThreshFilter::New();
            rthresh->SetInput(rlabels);
            rthresh->SetLowerThreshold(i);
            rthresh->SetUpperThreshold(i);
            rthresh->SetInsideValue(1);
            rthresh->SetOutsideValue(0);
            rthresh->Update();
            fname = prefix + suffix.str() + QI::OutExt();
            if (verbose)
                cout << "Writing output file " << fname << endl;
            QI::WriteImage(rimage, fname);
            fname = prefix + suffix.str() + "_mask" + QI::OutExt();
            if (verbose)
                cout << "Writing output mask " << fname << endl;
            QI::WriteImage(rthresh->GetOutput(), fname);
        }
    }
    return EXIT_SUCCESS;
}
