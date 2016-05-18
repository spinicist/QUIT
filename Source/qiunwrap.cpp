/*
 *  qiunwrap.cpp
 *
 *  Created by Tobias Wood on 11/06/2015.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <getopt.h>
#include <iostream>
#include "Eigen/Dense"

#include "itkImageSource.h"
#include "itkImageFileReader.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMultiplyImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkForwardFFTImageFilter.h"
#include "itkInverseFFTImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkFFTPadImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryCrossStructuringElement.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "QI/Types.h"
#include "QI/Util.h"

using namespace std;
using namespace Eigen;

namespace itk {

class DiscreteLaplacePhaseFilter : public ImageToImageFilter<QI::VolumeF, QI::VolumeF> {
protected:

public:
    /** Standard class typedefs. */
    typedef QI::VolumeF     TImage;

    typedef DiscreteLaplacePhaseFilter         Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef SmartPointer<Self>                 Pointer;
    typedef typename TImage::RegionType        RegionType;

    itkNewMacro(Self);
    itkTypeMacro(DiscreteLaplacePhaseFilter, DiscreteLaplacePhaseFilter);

    void SetInput(const TImage *img) override      { this->SetNthInput(0, const_cast<TImage*>(img)); }
    typename TImage::ConstPointer GetInput() const { return static_cast<const TImage *>(this->ProcessObject::GetInput(0)); }

    virtual void GenerateOutputInformation() override {
        Superclass::GenerateOutputInformation();
        auto op = this->GetOutput();
        op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
        op->Allocate();
    }

protected:
    DiscreteLaplacePhaseFilter() {
        this->SetNumberOfRequiredInputs(1);
        this->SetNumberOfRequiredOutputs(1);
        this->SetNthOutput(0, this->MakeOutput(0));
    }
    ~DiscreteLaplacePhaseFilter() {}

    DataObject::Pointer MakeOutput(unsigned int idx) {
        //std::cout <<  __PRETTY_FUNCTION__ << endl;
        if (idx == 0) {
            DataObject::Pointer output = (TImage::New()).GetPointer();
            return output.GetPointer();
        } else {
            std::cerr << "No output " << idx << std::endl;
            return NULL;
        }
    }

    virtual void ThreadedGenerateData(const RegionType &region, ThreadIdType threadId) override {
        //std::cout <<  __PRETTY_FUNCTION__ << endl;
        ConstNeighborhoodIterator<TImage>::RadiusType radius;
        radius.Fill(1);
        ConstNeighborhoodIterator<TImage> inputIter(radius, this->GetInput(), region);
        ImageRegionIterator<TImage> outputIter(this->GetOutput(), region);

        vector<ConstNeighborhoodIterator<TImage>::OffsetType> back, fwrd;
        back.push_back({{-1, 0, 0}});
        fwrd.push_back({{ 1, 0, 0}});
        back.push_back({{ 0,-1, 0}});
        fwrd.push_back({{ 0, 1, 0}});
        back.push_back({{ 0, 0,-1}});
        fwrd.push_back({{ 0, 0, 1}});

        TImage::SpacingType spacing = this->GetInput()->GetSpacing();
        TImage::SpacingType s_sqr = spacing * spacing;
        while(!inputIter.IsAtEnd()) {
            double sum = 0;
            double cphase = inputIter.GetCenterPixel();
            complex<double> c = std::polar(1., cphase);
            for (int i = 0; i < fwrd.size(); ++i) {
                double bphase = inputIter.GetPixel(back[i]);
                double fphase = inputIter.GetPixel(fwrd[i]);
                complex<double> b = std::polar(1., bphase);
                complex<double> f = std::polar(1., fphase);
                sum += std::arg((f*b)/(c*c)) / s_sqr[i];
            }
            outputIter.Set(sum / 7.);
            ++inputIter;
            ++outputIter;
        }
    }

private:
    DiscreteLaplacePhaseFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

class DiscreteInverseLaplace : public ImageSource<QI::VolumeF> {
public:
    typedef QI::VolumeF            TImage;
    typedef DiscreteInverseLaplace Self;
    typedef ImageSource<TImage>    Superclass;
    typedef SmartPointer<Self>     Pointer;

    itkNewMacro(Self);
    itkTypeMacro(DiscreteInverseLaplace, ImageSource);

    void SetImageProperties(const TImage *img) {
        m_region = img->GetLargestPossibleRegion();
        m_spacing = img->GetSpacing();
        m_direction = img->GetDirection();
        m_origin = img->GetOrigin();
    }

protected:
    typename TImage::RegionType    m_region;
    typename TImage::SpacingType   m_spacing;
    typename TImage::DirectionType m_direction;
    typename TImage::PointType     m_origin;

    DiscreteInverseLaplace(){}
    ~DiscreteInverseLaplace(){}
    virtual void GenerateData() override {
        typename TImage::Pointer output = this->GetOutput();
        output->SetRegions(m_region);
        output->Allocate();
        output->SetSpacing(m_spacing);
        output->SetDirection(m_direction);
        output->SetOrigin(m_origin);
        itk::ImageRegionIteratorWithIndex<TImage> imageIt(output,output->GetLargestPossibleRegion());
        imageIt.GoToBegin();
        imageIt.Set(0.); // There is a pole here
        ++imageIt;
        while(!imageIt.IsAtEnd()) {
            auto index = imageIt.GetIndex() - m_region.GetIndex(); // Might be padded to a negative start
            double val = 0;
            for (int i = 0; i < 3; i++) {
                val += 2. - 2. * cos(index[i] * 2. * M_PI / m_region.GetSize()[i]);
            }
            val /= 7.;
            imageIt.Set(1./val);
            ++imageIt;
        }
    }

private:
    DiscreteInverseLaplace(const Self &);
    void operator=(const Self &);
};

/*
 * Do it the hard way!
 */
class DiscreteInverseLaplace2 : public ImageSource<QI::VolumeF> {
public:
    typedef QI::VolumeF             TImage;
    typedef DiscreteInverseLaplace2 Self;
    typedef ImageSource<TImage>     Superclass;
    typedef SmartPointer<Self>      Pointer;

    itkNewMacro(Self);
    itkTypeMacro(DiscreteInverseLaplace2, ImageSource);

    void SetImageProperties(const TImage *img) {
        m_region = img->GetLargestPossibleRegion();
        m_spacing = img->GetSpacing();
        m_direction = img->GetDirection();
        m_origin = img->GetOrigin();
    }

protected:
    typename TImage::RegionType    m_region;
    typename TImage::SpacingType   m_spacing;
    typename TImage::DirectionType m_direction;
    typename TImage::PointType     m_origin;

    DiscreteInverseLaplace2(){}
    ~DiscreteInverseLaplace2(){}
    virtual void GenerateData() override {
        typename TImage::Pointer output = this->GetOutput();
        output->SetRegions(m_region);
        output->Allocate();
        output->SetSpacing(m_spacing);
        output->SetDirection(m_direction);
        output->SetOrigin(m_origin);

        typename TImage::Pointer filt = TImage::New();
        filt->SetRegions(m_region);
        filt->SetSpacing(m_spacing);
        filt->SetDirection(m_direction);
        filt->SetOrigin(m_origin);
        filt->Allocate();
        typename TImage::IndexType index = m_region.GetIndex();
        index[0] += m_region.GetSize()[0] / 2;
        index[1] += m_region.GetSize()[1] / 2;
        index[2] += m_region.GetSize()[2] / 2;
        filt->SetPixel(index, -6);
        index[0] -= 1; filt->SetPixel(index, 1);
        index[0] += 2; filt->SetPixel(index, 1); index[0] -= 1;
        index[1] -= 1; filt->SetPixel(index, 1);
        index[1] += 2; filt->SetPixel(index, 1); index[1] -= 1;
        index[2] -= 1; filt->SetPixel(index, 1);
        index[2] += 2; filt->SetPixel(index, 1); index[2] -= 1;
        auto shift = itk::FFTShiftImageFilter<TImage, TImage>::New();
        shift->SetInput(filt);
        typedef itk::ForwardFFTImageFilter<TImage> FFFTType;
        auto forwardFFT = FFFTType::New();
        forwardFFT->SetInput(shift->GetOutput());
        auto magFilt = itk::ComplexToModulusImageFilter<QI::VolumeXF, QI::VolumeF>::New();
        magFilt->SetInput(forwardFFT->GetOutput());
        magFilt->Update();
        auto divFilt = itk::DivideImageFilter<QI::VolumeF, QI::VolumeF, QI::VolumeF>::New();
        divFilt->SetConstant1(1.);
        divFilt->SetInput2(magFilt->GetOutput());
        divFilt->Update();
        QI::VolumeF::Pointer filt2 = divFilt->GetOutput();
        filt2->DisconnectPipeline();
        index = m_region.GetIndex();
        filt2->SetPixel(index,0);
        this->GraftOutput(filt2);
    }

private:
    DiscreteInverseLaplace2(const Self &);
    void operator=(const Self &);
};


} // End namespace itk

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: qiunwrap [options] input \n\
\n\
Input is a single wrapped phase volume\n\
\n\
Options:\n\
    --help, -h        : Print this message.\n\
    --verbose, -v     : Print more information.\n\
    --out, -o path    : Specify an output filename (default image base).\n\
    --mask, -m file   : Mask input with specified file.\n\
    --erode, -e N     : Erode mask by N mm (Default 1 mm).\n\
    --debug, -d       : Save all pipeline steps.\n\
    --threads, -T N   : Use N threads (default=hardware limit).\n"
};

const struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"verbose", no_argument, 0, 'v'},
    {"out", required_argument, 0, 'o'},
    {"mask", required_argument, 0, 'm'},
    {"erode", required_argument, 0, 'e'},
    {"debug", no_argument, 0, 'd'},
    {"threads", required_argument, 0, 'T'},
    {0, 0, 0, 0}
};
const char *short_options = "hvo:m:e:dT:";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();

    bool verbose = false, debug = false;
    float erodeRadius = 1;
    string prefix;
    QI::VolumeUC::Pointer mask = ITK_NULLPTR;
    int indexptr = 0, c;
    while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
        switch (c) {
            case 'v': verbose = true; break;
            case 'm': {
                if (verbose) cout << "Reading mask file " << optarg << endl;
                auto maskFile = itk::ImageFileReader<QI::VolumeF>::New();
                auto maskThresh = itk::BinaryThresholdImageFilter<QI::VolumeF, QI::VolumeUC>::New();
                maskFile->SetFileName(optarg);
                maskThresh->SetInput(maskFile->GetOutput());
                maskThresh->SetLowerThreshold(1.0);
                maskThresh->SetInsideValue(1);
                maskThresh->Update();
                mask = maskThresh->GetOutput();
                mask->DisconnectPipeline();
            } break;
            case 'e': erodeRadius = atof(optarg); break;
            case 'o':
                prefix = optarg;
                cout << "Output prefix will be: " << prefix << endl;
                break;
            case 'd': debug = true; break;
            case 'T': itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg)); break;
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
    if ((argc - optind) != 1) {
        cout << "Incorrect number of arguments." << endl << usage << endl;
        return EXIT_FAILURE;
    }
    if (verbose) cout << "Opening input file: " << argv[optind] << endl;
    string fname(argv[optind++]);
    if (prefix == "")
        prefix = fname.substr(0, fname.find(".nii"));
    string outname = prefix + "_unwrap" + QI::OutExt();
    if (verbose) cout << "Output filename: " << outname << endl;

    auto inFile = QI::ReadImage(fname);
    auto calcLaplace = itk::DiscreteLaplacePhaseFilter::New();
    calcLaplace->SetInput(inFile);
    calcLaplace->Update();
    if (debug) QI::WriteImage(calcLaplace->GetOutput(), prefix + "_step1_laplace" + QI::OutExt());

    QI::VolumeF::Pointer lap = calcLaplace->GetOutput();
    if (mask) {
        auto masker = itk::MaskImageFilter<QI::VolumeF, QI::VolumeUC>::New();
        masker->SetInput(calcLaplace->GetOutput());
        masker->SetMaskImage(mask);
        if (erodeRadius > 0) {
            typedef itk::BinaryBallStructuringElement<QI::VolumeUC::PixelType, 3> ElementType;
            ElementType structuringElement;
            ElementType::SizeType radii;
            auto spacing = mask->GetSpacing();
            radii[0] = ceil(erodeRadius / spacing[0]);
            radii[1] = ceil(erodeRadius / spacing[1]);
            radii[2] = ceil(erodeRadius / spacing[2]);
            structuringElement.SetRadius(radii);
            structuringElement.CreateStructuringElement();
            if (verbose) cout << "Eroding mask by " << erodeRadius << " mm (" << radii << " voxels)" << endl;
            typedef itk::BinaryErodeImageFilter <QI::VolumeUC, QI::VolumeUC, ElementType> BinaryErodeImageFilterType;
            BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
            erodeFilter->SetInput(mask);
            erodeFilter->SetErodeValue(1);
            erodeFilter->SetKernel(structuringElement);
            erodeFilter->Update();
            masker->SetMaskImage(erodeFilter->GetOutput());
            if (debug) QI::WriteImage(erodeFilter->GetOutput(), prefix + "_eroded_mask" + QI::OutExt());
        }
        if (verbose) cout << "Applying mask" << endl;
        masker->Update();
        lap = masker->GetOutput();
        lap->DisconnectPipeline();
        if (debug) QI::WriteImage(lap, prefix + "_step1_laplace_masked" + QI::OutExt());
    }

    if (verbose) cout << "Padding image to valid FFT size." << endl;
    typedef itk::FFTPadImageFilter<QI::VolumeF> PadFFTType;
    auto padFFT = PadFFTType::New();
    padFFT->SetInput(lap);
    padFFT->Update();
    if (debug) QI::WriteImage(padFFT->GetOutput(), prefix + "_step2_padFFT" + QI::OutExt());
    if (verbose) {
        cout << "Padded image size: " << padFFT->GetOutput()->GetLargestPossibleRegion().GetSize() << endl;
        cout << "Calculating Forward FFT." << endl;
    }
    typedef itk::ForwardFFTImageFilter<QI::VolumeF> FFFTType;
    auto forwardFFT = FFFTType::New();
    forwardFFT->SetInput(padFFT->GetOutput());
    forwardFFT->Update();
    if (debug) QI::WriteImage(forwardFFT->GetOutput(), prefix + "_step3_forwardFFT" + QI::OutExt());
    if (verbose) cout << "Generating Inverse Laplace Kernel." << endl;
    auto inverseLaplace = itk::DiscreteInverseLaplace::New();
    inverseLaplace->SetImageProperties(padFFT->GetOutput());
    inverseLaplace->Update();
    if (debug) QI::WriteImage(inverseLaplace->GetOutput(), prefix + "_inverse_laplace_filter" + QI::OutExt());
    if (verbose) cout << "Multiplying." << endl;
    auto mult = itk::MultiplyImageFilter<QI::VolumeXF, QI::VolumeF, QI::VolumeXF>::New();
    mult->SetInput1(forwardFFT->GetOutput());
    mult->SetInput2(inverseLaplace->GetOutput());
    if (debug) QI::WriteImage(mult->GetOutput(), prefix + "_step3_multFFT" + QI::OutExt());
    if (verbose) cout << "Inverse FFT." << endl;
    auto inverseFFT = itk::InverseFFTImageFilter<QI::VolumeXF, QI::VolumeF>::New();
    inverseFFT->SetInput(mult->GetOutput());
    inverseFFT->Update();
    if (debug) QI::WriteImage(inverseFFT->GetOutput(), prefix + "_step4_inverseFFT" + QI::OutExt());
    if (verbose) cout << "Extracting original size image" << endl;
    auto extract = itk::ExtractImageFilter<QI::VolumeF, QI::VolumeF>::New();
    extract->SetInput(inverseFFT->GetOutput());
    extract->SetDirectionCollapseToSubmatrix();
    extract->SetExtractionRegion(calcLaplace->GetOutput()->GetLargestPossibleRegion());
    extract->Update();
    if (debug) QI::WriteImage(extract->GetOutput(), prefix + "_step5_extract" + QI::OutExt());
    if (mask) {
        if (verbose) cout << "Re-applying mask" << endl;
        auto masker = itk::MaskImageFilter<QI::VolumeF, QI::VolumeUC>::New();
        masker->SetMaskImage(mask);
        masker->SetInput(extract->GetOutput());
        masker->Update();
        QI::WriteImage(masker->GetOutput(), outname);
    } else {
        QI::WriteImage(extract->GetOutput(), outname);
    }
    if (verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
