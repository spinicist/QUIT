/*
 *  qi_unwrap_laplace.cpp
 *
 *  Copyright (c) 2015, 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>

#include "itkImageSource.h"
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
#include "itkExtractImageFilter.h"

#include "ImageTypes.h"
#include "Util.h"
#include "ImageIO.h"
#include "Args.h"

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
    itkTypeMacro(Self, Superclass);

    void SetInput(const TImage *img) ITK_OVERRIDE      { this->SetNthInput(0, const_cast<TImage*>(img)); }
    typename TImage::ConstPointer GetInput() const { return static_cast<const TImage *>(this->ProcessObject::GetInput(0)); }

    void GenerateOutputInformation() ITK_OVERRIDE {
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

    DataObject::Pointer MakeOutput(ProcessObject::DataObjectPointerArraySizeType idx) override {
        //std::cout <<  __PRETTY_FUNCTION__ );
        if (idx == 0) {
            DataObject::Pointer output = (TImage::New()).GetPointer();
            return output.GetPointer();
        } else {
            std::cerr << "No output " << idx << std::endl;
            return nullptr;
        }
    }

    void ThreadedGenerateData(const RegionType &region, ThreadIdType /* Unused */) ITK_OVERRIDE {
        //std::cout <<  __PRETTY_FUNCTION__ );
        ConstNeighborhoodIterator<TImage>::RadiusType radius;
        radius.Fill(1);
        ConstNeighborhoodIterator<TImage> inputIter(radius, this->GetInput(), region);
        ImageRegionIterator<TImage> outputIter(this->GetOutput(), region);

        std::vector<ConstNeighborhoodIterator<TImage>::OffsetType> back, fwrd;
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
            std::complex<double> c = std::polar(1., cphase);
            for (size_t i = 0; i < fwrd.size(); ++i) {
                double bphase = inputIter.GetPixel(back[i]);
                double fphase = inputIter.GetPixel(fwrd[i]);
                std::complex<double> b = std::polar(1., bphase);
                std::complex<double> f = std::polar(1., fphase);
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
    void GenerateData() ITK_OVERRIDE {
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
    void GenerateData() ITK_OVERRIDE {
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
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Laplacian phase unwrapping\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "PHASE", "Wrapped phase image");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<int> erode(parser, "ERODE", "Erode mask by N mm (default 1)", {'e', "erode"}, 1);
    args::Flag debug(parser, "DEBUG", "Output debugging images", {'d', "debug"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI_LOG(verbose, "Opening input file: " << QI::CheckPos(input_path));
    auto inFile = QI::ReadImage(QI::CheckPos(input_path));
    std::string prefix = (outarg ? outarg.Get() : QI::StripExt(input_path.Get()));

    auto calcLaplace = itk::DiscreteLaplacePhaseFilter::New();
    calcLaplace->SetInput(inFile);
    calcLaplace->Update();
    if (debug) QI::WriteImage(calcLaplace->GetOutput(), prefix + "_step1_laplace" + QI::OutExt());

    QI::VolumeF::Pointer lap = calcLaplace->GetOutput();
    auto mask_img = mask ? QI::ReadImage<QI::VolumeUC>(mask.Get()) : ITK_NULLPTR;
    if (mask) {
        auto masker = itk::MaskImageFilter<QI::VolumeF, QI::VolumeUC>::New();
        masker->SetInput(calcLaplace->GetOutput());
        
        if (erode) {
            typedef itk::BinaryBallStructuringElement<QI::VolumeUC::PixelType, 3> ElementType;
            ElementType structuringElement;
            ElementType::SizeType radii;
            auto spacing = mask_img->GetSpacing();
            radii[0] = ceil(erode.Get() / spacing[0]);
            radii[1] = ceil(erode.Get() / spacing[1]);
            radii[2] = ceil(erode.Get() / spacing[2]);
            structuringElement.SetRadius(radii);
            structuringElement.CreateStructuringElement();
            QI_LOG(verbose, "Eroding mask by " << erode.Get() << " mm (" << radii << " voxels)" );
            typedef itk::BinaryErodeImageFilter <QI::VolumeUC, QI::VolumeUC, ElementType> BinaryErodeImageFilterType;
            BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
            erodeFilter->SetInput(mask_img);
            erodeFilter->SetErodeValue(1);
            erodeFilter->SetKernel(structuringElement);
            erodeFilter->Update();
            masker->SetMaskImage(erodeFilter->GetOutput());
            if (debug) QI::WriteImage(erodeFilter->GetOutput(), prefix + "_eroded_mask" + QI::OutExt());
        } else {
            masker->SetMaskImage(mask_img);
        }
        QI_LOG(verbose, "Applying mask" );
        masker->Update();
        lap = masker->GetOutput();
        lap->DisconnectPipeline();
        if (debug) QI::WriteImage(lap, prefix + "_step1_laplace_masked" + QI::OutExt());
    }

    QI_LOG(verbose, "Padding image to valid FFT size." );
    typedef itk::FFTPadImageFilter<QI::VolumeF> PadFFTType;
    auto padFFT = PadFFTType::New();
    padFFT->SetInput(lap);
    padFFT->Update();
    if (debug) QI::WriteImage(padFFT->GetOutput(), prefix + "_step2_padFFT" + QI::OutExt());
    QI_LOG(verbose, "Padded image size: " << padFFT->GetOutput()->GetLargestPossibleRegion().GetSize() <<
                    "Calculating Forward FFT.");
    typedef itk::ForwardFFTImageFilter<QI::VolumeF> FFFTType;
    auto forwardFFT = FFFTType::New();
    forwardFFT->SetInput(padFFT->GetOutput());
    forwardFFT->Update();
    if (debug) QI::WriteImage(forwardFFT->GetOutput(), prefix + "_step3_forwardFFT" + QI::OutExt());
    QI_LOG(verbose, "Generating Inverse Laplace Kernel.");
    auto inverseLaplace = itk::DiscreteInverseLaplace::New();
    inverseLaplace->SetImageProperties(padFFT->GetOutput());
    inverseLaplace->Update();
    if (debug) QI::WriteImage(inverseLaplace->GetOutput(), prefix + "_inverse_laplace_filter" + QI::OutExt());
    QI_LOG(verbose, "Multiplying." );
    auto mult = itk::MultiplyImageFilter<QI::VolumeXF, QI::VolumeF, QI::VolumeXF>::New();
    mult->SetInput1(forwardFFT->GetOutput());
    mult->SetInput2(inverseLaplace->GetOutput());
    if (debug) QI::WriteImage(mult->GetOutput(), prefix + "_step3_multFFT" + QI::OutExt());
    QI_LOG(verbose, "Inverse FFT." );
    auto inverseFFT = itk::InverseFFTImageFilter<QI::VolumeXF, QI::VolumeF>::New();
    inverseFFT->SetInput(mult->GetOutput());
    inverseFFT->Update();
    if (debug) QI::WriteImage(inverseFFT->GetOutput(), prefix + "_step4_inverseFFT" + QI::OutExt());
    QI_LOG(verbose, "Extracting original size image" );
    auto extract = itk::ExtractImageFilter<QI::VolumeF, QI::VolumeF>::New();
    extract->SetInput(inverseFFT->GetOutput());
    extract->SetDirectionCollapseToSubmatrix();
    extract->SetExtractionRegion(calcLaplace->GetOutput()->GetLargestPossibleRegion());
    extract->Update();
    if (debug) QI::WriteImage(extract->GetOutput(), prefix + "_step5_extract" + QI::OutExt());
    std::string outname = prefix + "_unwrap" + QI::OutExt();
    QI_LOG(verbose, "Output filename: " << outname );
    if (mask) {
        QI_LOG(verbose, "Re-applying mask" );
        auto masker = itk::MaskImageFilter<QI::VolumeF, QI::VolumeUC>::New();
        masker->SetMaskImage(mask_img);
        masker->SetInput(extract->GetOutput());
        masker->Update();
        QI::WriteImage(masker->GetOutput(), outname);
    } else {
        QI::WriteImage(extract->GetOutput(), outname);
    }
    QI_LOG(verbose, "Finished." );
    return EXIT_SUCCESS;
}
