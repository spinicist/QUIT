/*
 *  qicomplex.cpp
 *
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkMultiplyImageFilter.h"
#include "itkComplexToPhaseImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkComposeImageFilter.h"
#include "itkMagnitudeAndPhaseToComplexImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"

#include "Util.h"
#include "ImageIO.h"
#include "Args.h"

namespace itk {
template<typename TImage>
class NegateFilter : public InPlaceImageFilter<TImage> {
public:
    typedef NegateFilter                Self;
    typedef InPlaceImageFilter<TImage>  Superclass;
    typedef SmartPointer<Self>          Pointer;
    typedef typename TImage::RegionType        TRegion;
    typedef typename TImage::InternalPixelType TPixel;

    itkNewMacro(Self);
    itkTypeMacro(NegateFilter, InPlaceImageFilter);

private:
    NegateFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
    bool m_all = false, m_alternate = false;
protected:
    NegateFilter() {}
    ~NegateFilter() {}

public:
    void SetNegate(const bool all, const bool alternate) {
        m_all = all;
        m_alternate = alternate;
    }
    void ThreadedGenerateData(const TRegion &region, ThreadIdType threadId) ITK_OVERRIDE {
        typedef typename TImage::PixelType PixelType;
        ImageSliceConstIteratorWithIndex<TImage> inIt(this->GetInput(), region);
        ImageSliceIteratorWithIndex<TImage>      outIt(this->GetOutput(), region);

        inIt.SetFirstDirection(0);
        inIt.SetSecondDirection(1);
        outIt.SetFirstDirection(0);
        outIt.SetSecondDirection(1);
        inIt.GoToBegin();
        outIt.GoToBegin();
        PixelType mult(m_all ? -1 : 1, 0);
        while(!inIt.IsAtEnd()) {
            while (!inIt.IsAtEndOfSlice()) {
                while (!inIt.IsAtEndOfLine()) {
                    outIt.Set(mult * inIt.Get());
                    ++inIt;
                    ++outIt;
                }
                inIt.NextLine();
                outIt.NextLine();
            }
            if (m_alternate) {
                mult = mult * PixelType(-1, 0);
            }
            inIt.NextSlice();
            outIt.NextSlice();
        }
        // std::std::cout << "End " << __PRETTY_FUNCTION__ << std::std::endl;
    }
};

} // End namespace itk

/* Arguments defined here so they are available in the templated run function */
args::ArgumentParser parser(
    "Input is specified with lower case letters. A valid combination of inputs "
    " must be specified, e.g. real & imaginary or magnitude & phase.\n"
    "Output is specified with upper case letters. Any combination can be given\n"
    "http://github.com/spinicist/QUIT");

args::ValueFlag<std::string> in_mag(parser, "IN_MAG", "Input magnitude file", {'m', "mag"});
args::ValueFlag<std::string> in_pha(parser, "IN_PHA", "Input phase file", {'p', "pha"});
args::ValueFlag<std::string> in_real(parser, "IN_REAL", "Input real file", {'r', "real"});
args::ValueFlag<std::string> in_imag(parser, "IN_IMAG", "Input imaginary file", {'i', "imag"});
args::ValueFlag<std::string> in_complex(parser, "IN_CPLX", "Input complex file", {'x', "complex"});
args::ValueFlag<std::string> in_realimag(parser, "IN_REALIMAG", "Input real & imaginary file", {'l', "realimag"});

args::ValueFlag<std::string> out_mag(parser, "OUT_MAG", "Output magnitude file", {'M', "MAG"});
args::ValueFlag<std::string> out_pha(parser, "OUT_PHA", "Output phase file", {'P', "PHA"});
args::ValueFlag<std::string> out_real(parser, "OUT_REAL", "Output real file", {'R', "REAL"});
args::ValueFlag<std::string> out_imag(parser, "OUT_IMAG", "Output imaginary file", {'I', "IMAG"});
args::ValueFlag<std::string> out_complex(parser, "OUT_CPLX", "Output complex file", {'X', "COMPLEX"});

args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
args::Flag     use_double(parser, "DOUBLE", "Process & output at double precision", {'d', "double"});
args::Flag     fixge(parser, "FIX_GE", "Negate alternate slices (fixes lack of FFT shift)", {"fixge"});
args::Flag     negate(parser, "NEGATE", "Negate entire volume", {"negate"});

template<typename TPixel>
void Run() {
    typedef itk::Image<TPixel, 4>          TImage;
    typedef itk::Image<std::complex<TPixel>, 4> TXImage;
    typedef itk::ImageFileReader<TImage>   TReader;
    typedef itk::ImageFileReader<TXImage>  TXReader;
    typedef itk::ImageFileWriter<TImage>   TWriter;
    typedef itk::ImageFileWriter<TXImage>  TXWriter;

    typename TImage::Pointer img1 = ITK_NULLPTR, img2 = ITK_NULLPTR;
    typename TXImage::Pointer imgX = ITK_NULLPTR;

    if (in_real) {
        if (verbose) std::cout << "Reading real file: " << in_real.Get() << std::endl;
        img1 = QI::ReadImage<TImage>(in_real.Get());
        if (in_imag) {
            if (verbose) std::cout << "Reading imaginary file: " << in_imag.Get() << std::endl;
            img2 = QI::ReadImage<TImage>(in_imag.Get());
        } else {
            QI_FAIL("Must set real and imaginary inputs together");
        }
        auto compose = itk::ComposeImageFilter<TImage, TXImage>::New();
        compose->SetInput(0, img1);
        compose->SetInput(1, img2);
        compose->Update();
        imgX = compose->GetOutput();
        imgX->DisconnectPipeline();
    } else if (in_mag) {
        if (verbose) std::cout << "Reading magnitude file: " << in_mag.Get() << std::endl;
        img1 = QI::ReadImage<TImage>(in_mag.Get());
        if (in_pha) {
            if (verbose) std::cout << "Reading phase file: " << in_pha.Get() << std::endl;
            img2 = QI::ReadImage<TImage>(in_pha.Get());
        } else {
            QI_FAIL("Must set magnitude and phase inputs together");
        }
        auto compose = itk::MagnitudeAndPhaseToComplexImageFilter<TImage, TImage, TXImage>::New();
        compose->SetInput(0, img1);
        compose->SetInput(1, img2);
        compose->Update();
        imgX = compose->GetOutput();
        imgX->DisconnectPipeline();
    } else if (in_complex) {
        if (verbose) std::cout << "Reading complex file: " << in_complex.Get() << std::endl;
        imgX = QI::ReadImage<TXImage>(in_complex.Get());
    } else if (in_realimag) {
        if (verbose) std::cout << "Reading real/imaginary file: " << in_realimag.Get() << std::endl;
        auto img_both = QI::ReadImage<TImage>(in_realimag.Get());
        auto real_region = img_both->GetLargestPossibleRegion();
        auto imag_region = img_both->GetLargestPossibleRegion();
        real_region.GetModifiableSize()[3] = real_region.GetSize()[3] / 2;
        imag_region.GetModifiableSize()[3] = real_region.GetSize()[3];
        imag_region.GetModifiableIndex()[3] = real_region.GetSize()[3];
        auto extract_real = itk::RegionOfInterestImageFilter<TImage, TImage>::New();
        extract_real->SetRegionOfInterest(real_region);
        extract_real->SetInput(img_both);
        extract_real->Update();
        auto extract_imag = itk::RegionOfInterestImageFilter<TImage, TImage>::New();
        extract_imag->SetRegionOfInterest(imag_region);
        extract_imag->SetInput(img_both);
        extract_imag->Update();
        /* Hack round ITK changing the origin */
        extract_imag->GetOutput()->SetOrigin(extract_real->GetOutput()->GetOrigin());
        auto compose = itk::ComposeImageFilter<TImage, TXImage>::New();
        compose->SetInput(0, extract_real->GetOutput());
        compose->SetInput(1, extract_imag->GetOutput());
        compose->Update();
        imgX = compose->GetOutput();
        imgX->DisconnectPipeline();
    } else {
        QI_FAIL("No input files specified, use --help to see usage");
    }

    if (fixge || negate) {
        if (verbose && negate) std::cout << "Negating values" << std::endl;
        if (verbose && fixge)  std::cout << "Fixing GE lack of FFT-shift bug" << std::endl;
        auto fix = itk::NegateFilter<TXImage>::New();
        fix->SetInput(imgX);
        fix->SetNegate(negate, fixge);
        fix->Update();
        imgX = fix->GetOutput();
        imgX->DisconnectPipeline();
    }

    if (verbose) std::cout << "Writing output files" << std::endl;
    typename TWriter::Pointer write = TWriter::New();

    if (out_mag) {
        auto o = itk::ComplexToModulusImageFilter<TXImage, TImage>::New();
        o->SetInput(imgX);
        o->Update();
        QI::WriteImage(o->GetOutput(), out_mag.Get());
        if (verbose) std::cout << "Wrote magnitude image " << out_mag.Get() << std::endl;
    }
    if (out_pha) {
        auto o = itk::ComplexToPhaseImageFilter<TXImage, TImage>::New();
        o->SetInput(imgX);
        o->Update();
        QI::WriteImage(o->GetOutput(), out_pha.Get());
        if (verbose) std::cout << "Wrote phase image " << out_pha.Get() << std::endl;
    }
    if (out_real) {
        auto o = itk::ComplexToRealImageFilter<TXImage, TImage>::New();
        o->SetInput(imgX);
        o->Update();
        QI::WriteImage(o->GetOutput(), out_real.Get());
        if (verbose) std::cout << "Wrote real image " << out_real.Get() << std::endl;
    }
    if (out_imag) {
        auto o = itk::ComplexToImaginaryImageFilter<TXImage, TImage>::New();
        o->SetInput(imgX);
        o->Update();
        QI::WriteImage(o->GetOutput(), out_imag.Get());
        if (verbose) std::cout << "Wrote imaginary image " << out_imag.Get() << std::endl;
    }
    if (out_complex) {
        if (verbose) std::cout << "Writing complex image " << out_complex.Get() << std::endl;
        QI::WriteImage(imgX, out_complex.Get());
    }
}

int main(int argc, char **argv) {
    QI::ParseArgs(parser, argc, argv);
    if (use_double) {
        if (verbose) std::cout << "Using double precision" << std::endl;
        Run<double>();
    } else {
        if (verbose) std::cout << "Using float precision" << std::endl;
        Run<float>();
    }
    return EXIT_SUCCESS;
}
