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
#include <sstream>
#include "Eigen/Dense"

#include "itkImageSource.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMultiplyImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkComplexToComplexFFTImageFilter.h"
#include "itkInverseFFTImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkFFTPadImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkCastImageFilter.h"

#include "QI/Types.h"
#include "QI/Util.h"

using namespace std;
using namespace Eigen;

class FilterKernel {
    
};

class TukeyKernel : FilterKernel {
protected:
    double m_a = 0.75;
    double m_q = 0.25;
    int m_hx = 0, m_hy = 0, m_hz = 0;
    
public:
    void setSize(const int sx, const int sy, const int sz) {
        m_hx = sx / 2;  m_hy = sy / 2;  m_hz = sz / 2;
    }
    void setAQ(const double a, const double q) {
        m_a = a;
        m_q = q;
    }
    
    double operator() (const int x, const int y, const int z) const {
            const double rx = static_cast<double>(x - m_hx)/m_hx;
            const double ry = static_cast<double>(y - m_hy)/m_hy;
            const double rz = static_cast<double>(z - m_hz)/m_hz;
            const double r = sqrt((rx*rx + ry*ry + rz*rz) / 3);
            const double v = (r <= (1 - m_a)) ? 1 : 0.5*((1+m_q)+(1-m_q)*cos((M_PI/m_a)*(r - 1 + m_a)));
            return v;
    }
};


namespace itk {

template<typename ImageType>
class KSpaceFilter : public ImageToImageFilter<ImageType, ImageType> {
public:
    typedef ImageType                          TImage;
    typedef KSpaceFilter                       Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef SmartPointer<Self>                 Pointer;

protected:
    typedef typename TImage::RegionType::SizeType SizeType;
    typename TImage::RegionType    m_region;
    typename TImage::SpacingType   m_spacing;
    typename TImage::DirectionType m_direction;
    typename TImage::PointType     m_origin;

    TukeyKernel m_kernel;

    KSpaceFilter(){}
    ~KSpaceFilter(){}

public:
    itkNewMacro(Self);
    itkTypeMacro(KSpaceFilter, ImageSource);

    virtual void GenerateOutputInformation() override {
        Superclass::GenerateOutputInformation();
        SizeType size = this->GetOutput()->GetLargestPossibleRegion().GetSize();
        m_kernel.setSize(size[0], size[1], size[2]);
    }
    
    void SetKernel(const TukeyKernel &k) { m_kernel = k; }

protected:
    virtual void ThreadedGenerateData(const typename TImage::RegionType &region, ThreadIdType threadId) override {
        itk::ImageRegionIterator<TImage> outIt(this->GetOutput(),region);
        itk::ImageRegionConstIteratorWithIndex<TImage> inIt(this->GetInput(),region);
        inIt.GoToBegin();
        outIt.GoToBegin();
        while(!inIt.IsAtEnd()) {
            const auto index = inIt.GetIndex() - m_region.GetIndex(); // Might be padded to a negative start
            const typename TImage::PixelType::value_type k = m_kernel(index[0], index[1], index[2]);
            const typename TImage::PixelType v = inIt.Get();
            outIt.Set(v * k);
            ++inIt;
            ++outIt;
        }
    }

private:
    KSpaceFilter(const Self &);
    void operator=(const Self &);
};

} // End namespace itk

template<typename TImage> typename TImage::Pointer FilterVolume(const typename TImage::Pointer invol, const TukeyKernel &kernel, const bool verbose, const bool debug, const string &out_name);

template<> QI::ImageXD::Pointer FilterVolume<QI::ImageXD>(const typename QI::ImageXD::Pointer invol, const TukeyKernel &kernel, const bool verbose, const bool debug, const string &out_name) {
    if (verbose) cout << "Padding image to valid FFT size." << endl;
    typedef itk::FFTPadImageFilter<QI::ImageXD> TPad;
    typename TPad::Pointer padFFT = TPad::New();
    padFFT->SetInput(invol);
    padFFT->Update();
    if (debug) QI::WriteImage(padFFT->GetOutput(), out_name + "_step1_padFFT" + QI::OutExt());
    if (verbose) {
        cout << "Padded image size: " << padFFT->GetOutput()->GetLargestPossibleRegion().GetSize() << endl;
        cout << "Calculating Forward FFT." << endl;
    }
    
    typedef itk::ComplexToComplexFFTImageFilter<QI::ImageXD> FFTType;
    auto forwardFFT = FFTType::New();
    forwardFFT->SetInput(padFFT->GetOutput());
    forwardFFT->Update();

    typedef itk::FFTShiftImageFilter<QI::ImageXD, QI::ImageXD> TShift;
    TShift::Pointer shiftFFT = TShift::New();
    shiftFFT->SetInput(forwardFFT->GetOutput());
    shiftFFT->Update();
    if (debug) QI::WriteImage(shiftFFT->GetOutput(), out_name + "_step2_forwardFFT" + QI::OutExt());

    if (verbose) cout << "Filtering." << endl;
    typedef itk::KSpaceFilter<QI::ImageXD> TFilter;
    TFilter::Pointer k_filter = TFilter::New();
    k_filter->SetKernel(kernel);
    k_filter->SetInput(shiftFFT->GetOutput());
    if (debug) QI::WriteImage   (k_filter->GetOutput(), out_name + "_step3_multFFT" + QI::OutExt());

    auto unshiftFFT = TShift::New();
    unshiftFFT->SetInput(k_filter->GetOutput());
    unshiftFFT->Update();

    if (verbose) cout << "Inverse FFT." << endl;
    auto inverseFFT = FFTType::New();
    inverseFFT->SetTransformDirection(FFTType::INVERSE);
    inverseFFT->SetInput(unshiftFFT->GetOutput());
    inverseFFT->Update();

    if (debug) QI::WriteImage(inverseFFT->GetOutput(), out_name + "_step4_inverseFFT" + QI::OutExt());
    if (verbose) cout << "Extracting original size image" << endl;
    auto extract = itk::ExtractImageFilter<QI::ImageXD, QI::ImageXD>::New();
    extract->SetInput(inverseFFT->GetOutput());
    extract->SetDirectionCollapseToSubmatrix();
    extract->SetExtractionRegion(invol->GetLargestPossibleRegion());
    extract->Update();
    typename QI::ImageXD::Pointer outvol = extract->GetOutput();
    outvol->DisconnectPipeline();
    return outvol;
}

template<> QI::ImageD::Pointer FilterVolume<QI::ImageD>(const typename QI::ImageD::Pointer invol, const TukeyKernel &kernel, const bool verbose, const bool debug, const string &out_name) {
    typedef itk::CastImageFilter<QI::ImageD, QI::ImageXD> TCast;
    auto caster = TCast::New();
    caster->SetInput(invol);
    caster->Update();
    QI::ImageXD::Pointer filtered = FilterVolume<QI::ImageXD>(caster->GetOutput(), kernel, verbose, debug, out_name);
    if (verbose) cout << "Taking complex modulus" << endl;
    typedef itk::ComplexToModulusImageFilter<QI::ImageXD, QI::ImageD> TMod;
    auto mod = TMod::New();
    mod->SetInput(filtered);
    mod->Update();
    QI::ImageD::Pointer final = mod->GetOutput();
    final->DisconnectPipeline();
    return final;
}


template<typename TPixel>
void run_pipeline(const std::string &in_name, const std::string &out_name, const bool verbose, const bool debug, const TukeyKernel &kernel) {
    typedef itk::Image<TPixel, 4> TSeries;
    typedef itk::Image<TPixel, 3> TVolume;

    if (verbose) cout << "Opening input file: " << in_name << endl;
    typename TSeries::Pointer vols = QI::ReadImage<TSeries>(in_name);
    auto region = vols->GetLargestPossibleRegion();
    const size_t nvols = region.GetSize()[3]; // Save for the loop
    region.GetModifiableSize()[3] = 0;
    auto extract = itk::ExtractImageFilter<TSeries, TVolume>::New();
    extract->SetInput(vols);
    extract->SetDirectionCollapseToSubmatrix();
    typename itk::PasteImageFilter<TSeries>::Pointer paster = itk::PasteImageFilter<TSeries>::New();
    typename itk::CastImageFilter<TVolume, TSeries>::Pointer caster = itk::CastImageFilter<TVolume, TSeries>::New();
    paster->SetDestinationImage(vols);
    //paster->InPlaceOn();
    for (int i = 0; i < nvols; i++) {
        region.GetModifiableIndex()[3] = i;
        extract->SetExtractionRegion(region);
        if (verbose) cout << "Extracting volume " << i << endl;
        cout << extract->GetInPlace() << endl;
        extract->Update();
        if (verbose) cout << "Launching filter pipeline" << endl;
        typename TVolume::Pointer outvol = FilterVolume<TVolume>(extract->GetOutput(), kernel, verbose, debug, out_name);
        if (verbose) cout << "Casting to 4D" << endl;
        caster->SetInput(outvol);
        caster->Update();
        if (verbose) cout << "Pasting" << endl;
        paster->SetSourceImage(caster->GetOutput());
        paster->SetSourceRegion(caster->GetOutput()->GetLargestPossibleRegion());
        paster->SetDestinationImage(vols);
        paster->SetDestinationIndex(region.GetIndex());
        paster->Update();
        vols = paster->GetOutput();
    }
    cout << "paster in place " << paster->GetInPlace() << endl;
    if (verbose) cout << "Writing output:" << out_name + QI::OutExt() << endl;
    QI::WriteImage(vols, out_name + QI::OutExt());
    if (verbose) cout << "Finished." << endl;
}

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: qikfilter [options] input \n\
\n\
Filter images in k-Space\n\
\n\
Options:\n\
    --help, -h           : Print this message.\n\
    --verbose, -v        : Print more information.\n\
    --out, -o path       : Specify an output filename (default image base).\n\
    --filter, -f \"a q\" : Set Filter a and q values (default 0.75 0.25).\n\
    --debug, -d          : Save all pipeline steps.\n\
    --complex, -x        : Input/output data is complex.\n\
    --threads, -T N      : Use N threads (default=hardware limit).\n"
};

const struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"verbose", no_argument, 0, 'v'},
    {"out", required_argument, 0, 'o'},
    {"filter", required_argument, 0, 'f'},
    {"debug", no_argument, 0, 'd'},
    {"complex", no_argument, 0, 'x'},
    {"threads", required_argument, 0, 'T'},
    {0, 0, 0, 0}
};
const char *short_options = "hvo:m:e:dxT:";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();

    bool verbose = false, debug = false, is_complex = false;
    double filter_a = 0.75, filter_q = 0.25;
    string out_name;
    int indexptr = 0, c;
    while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
        switch (c) {
            case 'v': verbose = true; break;
            case 'o':
                out_name = optarg;
                cout << "Output filename will be: " << out_name << endl;
                break;
            case 'f': {
                stringstream f(optarg);
                f >> filter_a >> filter_q;
            } break;
            case 'd': debug = true; break;
            case 'x': is_complex = true; break;
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

    string in_name(argv[optind++]);
    if (out_name == "")
        out_name = in_name.substr(0, in_name.find(".")) + "_filtered";

    TukeyKernel filter_kernel;
    filter_kernel.setAQ(filter_a, filter_q);
    if (is_complex) {
        if (verbose) cout << "Data is complex." << endl;
        run_pipeline<complex<double>>(in_name, out_name, verbose, debug, filter_kernel);
    } else {
        if (verbose) cout << "Data is real." << endl;
        run_pipeline<double>(in_name, out_name, verbose, debug, filter_kernel);
    }

    return EXIT_SUCCESS;
}
