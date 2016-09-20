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

#include <getopt.h>
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

#include "QI/Util.h"

using namespace std;

namespace itk {
template<typename TImage>
class FixGEFilter : public InPlaceImageFilter<TImage> {
public:
    typedef FixGEFilter                 Self;
    typedef InPlaceImageFilter<TImage>  Superclass;
    typedef SmartPointer<Self>          Pointer;
    typedef typename TImage::RegionType        TRegion;
    typedef typename TImage::InternalPixelType TPixel;

    itkNewMacro(Self);
    itkTypeMacro(FixGEFilter, InPlaceImageFilter);

private:
    FixGEFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented

protected:
    FixGEFilter() {}
    ~FixGEFilter() {}

public:
    virtual void ThreadedGenerateData(const TRegion &region, ThreadIdType threadId) ITK_OVERRIDE {
        //std::cout <<  __PRETTY_FUNCTION__ << std::endl;
        typedef typename TImage::PixelType PixelType;
        ImageSliceConstIteratorWithIndex<TImage> inIt(this->GetInput(), region);
        ImageSliceIteratorWithIndex<TImage>      outIt(this->GetOutput(), region);

        inIt.SetFirstDirection(0);
        inIt.SetSecondDirection(1);
        outIt.SetFirstDirection(0);
        outIt.SetSecondDirection(1);
        inIt.GoToBegin();
        outIt.GoToBegin();
        PixelType mult(1, 0);
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
            mult = mult * PixelType(-1, 0);
            inIt.NextSlice();
            outIt.NextSlice();
        }
        // std::cout << "End " << __PRETTY_FUNCTION__ << std::endl;
    }
};

} // End namespace itk

const string usage {
"Usage is: qcomplex [input options] [output options] [other options] \n\
\n\
Input is specified with lower case letters. One of the following\n\
combinations must be specified:\n\
    -m mag_image -p phase_image\n\
    -r real_image -i imaginary_image\n\
    -x complex_image\n\
\n\
Output is specified with upper case letters. One or more of the\n\
following can be specified:\n\
    -M : Output a magnitude image\n\
    -P : Output a phase image\n\
    -R : Output a real image\n\
    -I : Output an imaginary image\n\
    -X : Output a complex image\n\
\n\
Other options:\n\
    --negate : Multiply everything by -1 before output.\n\
    --double : Use double precision instead of float\n\
    --fixge  : Fix alternate slice problem with GE data.\n\
\n\
Example:\n\
    qicomplex -m mag.nii -p phase.nii -R real.nii -I imag.nii\n"
};
int verbose = false, use_double = false, fixge = false, do_negate = false;
const struct option long_options[] =
{
    {"help", no_argument, 0, 'h'},
    {"verbose", no_argument, 0, 'v'},
    {"double", no_argument, &use_double, 1},
    {"fixge", no_argument, &fixge, 1},
    {"negate", no_argument, &do_negate, 1},
    {0, 0, 0, 0}
};
const char* short_options = "hvm:M:p:P:r:R:i:I:x:X:";
int c, index_ptr = 0;

template<typename TPixel> void Run(int argc, char **argv) {
    typedef itk::Image<TPixel, 4>          TImage;
    typedef itk::Image<complex<TPixel>, 4> TXImage;
    typedef itk::ImageFileReader<TImage>   TReader;
    typedef itk::ImageFileReader<TXImage>  TXReader;
    typedef itk::ImageFileWriter<TImage>   TWriter;
    typedef itk::ImageFileWriter<TXImage>  TXWriter;

    typename TImage::Pointer img1 = ITK_NULLPTR, img2 = ITK_NULLPTR;
    typename TXImage::Pointer imgX = ITK_NULLPTR;

    if (verbose) cout << "Reading input files" << endl;
    bool ri = false, x = false;
    optind = 1;
    while ((c = getopt_long(argc, argv, short_options,  long_options, &index_ptr)) != -1) {
        switch (c) {
            case 'r':
                ri = true;
            case 'm': {
                auto read = TReader::New();
                read->SetFileName(optarg);
                read->Update();
                img1 = read->GetOutput();
            } break;
            case 'i':
                ri = true;
            case 'p': {
                auto read = TReader::New();
                read->SetFileName(optarg);
                read->Update();
                img2 = read->GetOutput();
            } break;
            case 'x': {
                x = true;
                typename TXReader::Pointer readX = TXReader::New();
                readX->SetFileName(optarg);
                readX->Update();
                imgX = readX->GetOutput();
            } break;
            default: break;
        }
    }

    if (x) {
        // Nothing to see here
    } else if (ri) {
        if (!(img1 && img2)) {
            QI_EXCEPTION("Must set real and imaginary inputs");
        }
        if (verbose) cout << "Combining real and imaginary input" << endl;
        auto compose = itk::ComposeImageFilter<TImage, TXImage>::New();
        compose->SetInput(0, img1);
        compose->SetInput(1, img2);
        compose->Update();
        imgX = compose->GetOutput();
    } else {
        if (!(img1 && img2)) {
            QI_EXCEPTION("Must set magnitude and phase inputs");
        }
        if (verbose) cout << "Combining magnitude and phase input" << endl;
        auto compose = itk::MagnitudeAndPhaseToComplexImageFilter<TImage, TImage, TXImage>::New();
        compose->SetInput(0, img1);
        compose->SetInput(1, img2);
        compose->Update();
        imgX = compose->GetOutput();
    }

    if (fixge) {
        if (verbose) cout << "Fixing GE phase bug" << endl;
        auto fix = itk::FixGEFilter<TXImage>::New();
        fix->SetInput(imgX);
        fix->Update();
        imgX = fix->GetOutput();
    }

    if (do_negate) {
        if (verbose) cout << "Negating values" << endl;
        auto neg = itk::MultiplyImageFilter<TXImage, TXImage, TXImage>::New();
        neg->SetInput(imgX);
        neg->SetConstant(complex<TPixel>(-1.0, 0.0));
        neg->Update();
        imgX = neg->GetOutput();
    }
    if (verbose) cout << "Writing output files" << endl;
    typename TWriter::Pointer write = TWriter::New();
    optind = 1;
    while ((c = getopt_long(argc, argv, short_options, long_options, &index_ptr)) != -1) {
        switch (c) {
            case 'M': {
                auto o = itk::ComplexToModulusImageFilter<TXImage, TImage>::New();
                o->SetInput(imgX);
                write->SetFileName(optarg);
                write->SetInput(o->GetOutput());
                write->Update();
                if (verbose) cout << "Wrote magnitude image " + string(optarg) << endl;
            } break;
            case 'P': {
                auto o = itk::ComplexToPhaseImageFilter<TXImage, TImage>::New();
                o->SetInput(imgX);
                write->SetFileName(optarg);
                write->SetInput(o->GetOutput());
                write->Update();
                if (verbose) cout << "Wrote phase image " + string(optarg) << endl;
            } break;
            case 'R': {
                auto o = itk::ComplexToRealImageFilter<TXImage, TImage>::New();
                o->SetInput(imgX);
                write->SetFileName(optarg);
                write->SetInput(o->GetOutput());
                write->Update();
                if (verbose) cout << "Wrote real image " + string(optarg) << endl;
            } break;
            case 'I': {
                auto o = itk::ComplexToImaginaryImageFilter<TXImage, TImage>::New();
                o->SetInput(imgX);
                write->SetFileName(optarg);
                write->SetInput(o->GetOutput());
                write->Update();
                if (verbose) cout << "Wrote imaginary image " + string(optarg) << endl;
            } break;
            case 'X': {
                auto writeX = TXWriter::New();
                writeX->SetFileName(optarg);
                writeX->SetInput(imgX);
                writeX->Update();
                if (verbose) cout << "Wrote complex image " + string(optarg) << endl;
            } break;
            default: break;
        }
    }
}

int main(int argc, char **argv) {
    // Do one pass for the general options
    bool have_some_options = false;
    while ((c = getopt_long(argc, argv, short_options, long_options, &index_ptr)) != -1) {
        switch (c) {
            case 'v': verbose = true; break;
            case 'h':
                cout << QI::GetVersion() << endl << usage << endl;
                return EXIT_SUCCESS;
            case 0: break; // A flag
            case 'm': case 'M': case 'p': case 'P':
            case 'r': case 'R': case 'i': case 'I':
            case 'x': case 'X':
                have_some_options = true;
                break;
            case '?': // getopt will print an error message
                return EXIT_FAILURE;
            default:
            cout << "Unhandled option " << string(1, c) << endl;
            return EXIT_FAILURE;
        }
    }

    if (!have_some_options) {
        cout << QI::GetVersion() << endl << usage << endl;
        return EXIT_FAILURE;
    }

    if (use_double) {
        if (verbose) cout << "Using double precision" << endl;
        Run<double>(argc, argv);
    } else {
        if (verbose) cout << "Using float precision" << endl;
        Run<float>(argc, argv);
    }

    return EXIT_SUCCESS;
}
