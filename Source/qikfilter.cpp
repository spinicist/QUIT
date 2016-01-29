/*
 *  qikfilter.cpp
 *
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <time.h>
#include <getopt.h>
#include <iostream>
#include <atomic>
#include "Eigen/Dense"

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

#include "Types.h"
#include "Util.h"

using namespace std;
using namespace Eigen;

namespace itk {

class TukeyFilter : public ImageSource<QI::ImageF> {
public:
    typedef QI::ImageF          TImage;
    typedef TukeyFilter         Self;
    typedef ImageSource<TImage> Superclass;
    typedef SmartPointer<Self>  Pointer;

    itkNewMacro(Self);
    itkTypeMacro(TukeyFilter, ImageSource);

    void SetImageProperties(const TImage *img) {
        m_region = img->GetLargestPossibleRegion();
        m_spacing = img->GetSpacing();
        m_direction = img->GetDirection();
        m_origin = img->GetOrigin();
    }

protected:
    typedef typename TImage::RegionType::SizeType SizeType;
    typename TImage::RegionType    m_region;
    typename TImage::SpacingType   m_spacing;
    typename TImage::DirectionType m_direction;
    typename TImage::PointType     m_origin;

    double m_a = 0.75;
    double m_q = 0.25;

    TukeyFilter(){}
    ~TukeyFilter(){}
    virtual void GenerateData() override {
        typename TImage::Pointer output = this->GetOutput();
        output->SetRegions(m_region);
        output->Allocate();
        output->SetSpacing(m_spacing);
        output->SetDirection(m_direction);
        output->SetOrigin(m_origin);

        itk::ImageRegionIteratorWithIndex<TImage> imageIt(output,output->GetLargestPossibleRegion());
        SizeType size = m_region.GetSize();
        int hx = size[0] / 2; int hy = size[1] / 2; int hz = size[2] / 2;

        const double rad_k = sqrt(static_cast<double>((hx*hx)+(hy*hy)+(hz*hz)));
        imageIt.GoToBegin();
        ++imageIt;
        while(!imageIt.IsAtEnd()) {
            const auto index = imageIt.GetIndex() - m_region.GetIndex(); // Might be padded to a negative start
            int x = index[0]; int y = index[1]; int z = index[2];
            const double r = sqrt(static_cast<double>((x-hx)*(x-hx)+(y-hy)*(y-hy)+(z-hz)*(z-hz))) / rad_k;
            const double v = (r <= (1 - m_a)) ? 1 : 0.5*((1+m_q)+(1-m_q)*cos((M_PI/m_a)*(r - 1 + m_a)));
            imageIt.Set(v);
            ++imageIt;
        }
    }

private:
    TukeyFilter(const Self &);
    void operator=(const Self &);
};

} // End namespace itk

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: qikfilter [options] input \n\
\n\
Filter images in k-Space\n\
\n\
Options:\n\
    --help, -h        : Print this message.\n\
    --verbose, -v     : Print more information.\n\
    --out, -o path    : Specify an output filename (default image base).\n\
    --debug, -d       : Save all pipeline steps.\n\
    --threads, -T N   : Use N threads (default=hardware limit).\n"
};

const struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"verbose", no_argument, 0, 'v'},
    {"out", required_argument, 0, 'o'},
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
    string prefix;
    int indexptr = 0, c;
    while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
        switch (c) {
            case 'v': verbose = true; break;
            case 'o':
                prefix = optarg;
                cout << "Output prefix will be: " << prefix << endl;
                break;
            case 'd': debug = true; break;
            case 'T': itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg)); break;
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
    if ((argc - optind) != 1) {
        cout << "Incorrect number of arguments." << endl << usage << endl;
        return EXIT_FAILURE;
    }
    if (verbose) cout << "Opening input file: " << argv[optind] << endl;
    string fname(argv[optind++]);
    if (prefix == "")
        prefix = fname.substr(0, fname.find(".nii"));
    string outname = prefix + "_filtered" + QI::OutExt();
    if (verbose) cout << "Output filename: " << outname << endl;

    auto inFile = QI::ImageReaderF::New();
    inFile->SetFileName(fname);
    inFile->Update(); // Need the size info

    if (verbose) cout << "Padding image to valid FFT size." << endl;
    typedef itk::FFTPadImageFilter<QI::ImageF> PadFFTType;
    auto padFFT = PadFFTType::New();
    padFFT->SetInput(inFile->GetOutput());
    padFFT->Update();
    if (debug) QI::WriteImage(padFFT->GetOutput(), prefix + "_step1_padFFT" + QI::OutExt());
    if (verbose) {
        cout << "Padded image size: " << padFFT->GetOutput()->GetLargestPossibleRegion().GetSize() << endl;
        cout << "Calculating Forward FFT." << endl;
    }
    typedef itk::ForwardFFTImageFilter<QI::ImageF> FFFTType;
    auto forwardFFT = FFFTType::New();
    forwardFFT->SetInput(padFFT->GetOutput());
    forwardFFT->Update();
    if (debug) QI::WriteImage<QI::ImageXF>(forwardFFT->GetOutput(), prefix + "_step2_forwardFFT" + QI::OutExt());

    typedef itk::FFTShiftImageFilter<QI::ImageXF, QI::ImageXF> ShiftType;
    auto shiftFFT = ShiftType::New();
    shiftFFT->SetInput(forwardFFT->GetOutput());
    shiftFFT->Update();

    if (verbose) cout << "Generating K-Space Filter." << endl;
    auto k_filter = itk::TukeyFilter::New();
    k_filter->SetImageProperties(padFFT->GetOutput());
    k_filter->Update();
    if (debug) QI::WriteImage(k_filter->GetOutput(), prefix + "_kfilter" + QI::OutExt());
    if (verbose) cout << "Multiplying." << endl;
    auto mult = itk::MultiplyImageFilter<QI::ImageXF, QI::ImageF, QI::ImageXF>::New();
    mult->SetInput1(shiftFFT->GetOutput());
    mult->SetInput2(k_filter->GetOutput());
    if (debug) QI::WriteImage<QI::ImageXF>(mult->GetOutput(), prefix + "_step3_multFFT" + QI::OutExt());

    auto unshiftFFT = ShiftType::New();
    unshiftFFT->SetInput(mult->GetOutput());
    unshiftFFT->Update();

    if (verbose) cout << "Inverse FFT." << endl;
    auto inverseFFT = itk::InverseFFTImageFilter<QI::ImageXF, QI::ImageF>::New();
    inverseFFT->SetInput(unshiftFFT->GetOutput());
    inverseFFT->Update();
    if (debug) QI::WriteImage(inverseFFT->GetOutput(), prefix + "_step4_inverseFFT" + QI::OutExt());
    if (verbose) cout << "Extracting original size image" << endl;
    auto extract = itk::ExtractImageFilter<QI::ImageF, QI::ImageF>::New();
    extract->SetInput(inverseFFT->GetOutput());
    extract->SetDirectionCollapseToSubmatrix();
    extract->SetExtractionRegion(inFile->GetOutput()->GetLargestPossibleRegion());
    extract->Update();
    auto outFile = QI::ImageWriterF::New();
    outFile->SetInput(extract->GetOutput());
    outFile->SetFileName(outname);
    if (verbose) cout << "Writing output." << endl;
    outFile->Update();
    if (verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
