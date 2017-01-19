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

#include <memory>
#include <getopt.h>
#include <iostream>
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
#include "QI/Kernels.h"
#include "QI/IO.h"

using namespace std;
using namespace Eigen;

namespace itk {

template<typename ImageType>
class KSpaceFilter : public ImageToImageFilter<ImageType, ImageType> {
public:
    typedef ImageType                          TImage;
    typedef KSpaceFilter                       Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef SmartPointer<Self>                 Pointer;

protected:
    shared_ptr<QI::FilterKernel> m_kernel;
    bool m_WriteKernel = false;

    KSpaceFilter(){}
    ~KSpaceFilter(){}

public:
    itkNewMacro(Self);
    itkTypeMacro(KSpaceFilter, ImageSource);
    itkSetMacro(WriteKernel, bool);

    void SetKernel(const shared_ptr<QI::FilterKernel> &k) { m_kernel = k; }
protected:
    virtual void ThreadedGenerateData(const typename TImage::RegionType &region, ThreadIdType threadId) ITK_OVERRIDE {
        typename TImage::RegionType fullRegion = this->GetInput()->GetLargestPossibleRegion();
        typename TImage::IndexType startIndex = fullRegion.GetIndex();
        typename TImage::SizeType S = fullRegion.GetSize();
        itk::ImageRegionIterator<TImage> outIt(this->GetOutput(),region);
        itk::ImageRegionConstIteratorWithIndex<TImage> inIt(this->GetInput(),region);
        inIt.GoToBegin();
        outIt.GoToBegin();
        while(!inIt.IsAtEnd()) {
            const auto I = inIt.GetIndex() - startIndex; // Might be padded to a negative start
            const double rx = fmod(static_cast<double>(2.*I[0])/S[0] + 1.0, 2.0) - 1.0;
            const double ry = fmod(static_cast<double>(2.*I[1])/S[1] + 1.0, 2.0) - 1.0;
            const double rz = fmod(static_cast<double>(2.*I[2])/S[2] + 1.0, 2.0) - 1.0;
            const double r = sqrt((rx*rx + ry*ry + rz*rz) / 3);
            const typename TImage::PixelType::value_type k = m_kernel->value(r);
            const typename TImage::PixelType v = inIt.Get();
            if (m_WriteKernel) {
                outIt.Set(k);
            } else {
                outIt.Set(v * k);
            }
            ++inIt;
            ++outIt;
        }
    }

private:
    KSpaceFilter(const Self &);
    void operator=(const Self &);
};

} // End namespace itk

QI::SeriesXF::Pointer run_pipeline(QI::SeriesXF::Pointer vols, const bool verbose, const int debug, const shared_ptr<QI::FilterKernel> &kernel) {
    typedef itk::ExtractImageFilter<QI::SeriesXF, QI::VolumeXD> TExtract;
    typedef itk::TileImageFilter<QI::VolumeXF, QI::SeriesXF>    TTile;
    typedef itk::CastImageFilter<QI::VolumeXD, QI::VolumeXF>    TCast;
    typedef itk::FFTPadImageFilter<QI::VolumeXD>                TPad;
    typedef itk::ComplexToComplexFFTImageFilter<QI::VolumeXD>   TFFT;
    typedef itk::KSpaceFilter<QI::VolumeXD>                     TFilter;
    typedef itk::ExtractImageFilter<QI::VolumeXD, QI::VolumeXF> TUnpad;
    
    auto region = vols->GetLargestPossibleRegion();
    const size_t nvols = region.GetSize()[3]; // Save for the loop
    region.GetModifiableSize()[3] = 0;
    QI::VolumeXD::RegionType unpad_region;
    for (int i = 0; i < 3; i++) {
        unpad_region.GetModifiableSize()[i] = region.GetSize()[i];
    }

    itk::FixedArray<unsigned int, 4> layout;
    layout[0] = layout[1] = layout[2] = 1;
    layout[3] = nvols;

    auto extract = TExtract::New();
    auto pad     = TPad::New();
    auto forward = TFFT::New();
    auto k_filter = TFilter::New();
    // inverse is declared in the loop due to a weird bug
    auto unpadder = TUnpad::New();
    auto tile   = TTile::New();
    
    extract->SetInput(vols);
    extract->SetDirectionCollapseToSubmatrix();
    pad->SetInput(extract->GetOutput());
    forward->SetInput(pad->GetOutput());
    k_filter->SetInput(forward->GetOutput());
    k_filter->SetKernel(kernel);
    unpadder->SetDirectionCollapseToSubmatrix();
    unpadder->SetExtractionRegion(unpad_region);
    tile->SetLayout(layout);
    
    for (int i = 0; i < nvols; i++) {
        region.GetModifiableIndex()[3] = i;
        if (verbose) cout << "Processing volume " << i << endl;
        
        extract->SetExtractionRegion(region);
        extract->Update();
        auto inverse = TFFT::New();
        inverse->SetTransformDirection(TFFT::INVERSE);
        inverse->SetInput(k_filter->GetOutput());
        unpadder->SetInput(inverse->GetOutput());
        inverse->Update(); // Need a separate update to avoid region bug
        unpadder->Update();
        QI::VolumeXF::Pointer v = unpadder->GetOutput();
        tile->SetInput(i, v);
        v->DisconnectPipeline();
        if (debug) {
            cout << "Writing debug images for volume " << i << endl;
            QI::WriteMagnitudeImage(extract->GetOutput(), "kextract_" + to_string(i) + ".nii");
            QI::WriteMagnitudeImage(pad->GetOutput(), "kpad_" + to_string(i) + ".nii");
            QI::WriteMagnitudeImage(forward->GetOutput(), "kforward_" + to_string(i) + ".nii");
            QI::WriteMagnitudeImage(k_filter->GetOutput(), "kfiltered_" + to_string(i) + ".nii");
            k_filter->SetWriteKernel(true);
            k_filter->Update();
            QI::WriteMagnitudeImage(k_filter->GetOutput(), "kfilter_" + to_string(i) + ".nii");
            QI::WriteMagnitudeImage(inverse->GetOutput(), "kinverse_" + to_string(i) + ".nii");
        }
    }
    if (verbose) cout << "Finished." << endl;
    tile->Update();
    auto dir = vols->GetDirection();
    auto spc = vols->GetSpacing();
    vols = tile->GetOutput();
    vols->SetDirection(dir);
    vols->SetSpacing(spc);
    vols->DisconnectPipeline();
    return vols;
}

//******************************************************************************
// Main
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
    --filter, -f F       : Choose a filter (see below).\n\
    --debug, -d          : Save all pipeline steps.\n\
    --complex_in         : Read complex data.\n\
    --complex_out        : Output complex data.\n\
    --threads, -T N      : Use N threads (default=hardware limit).\n\
\n\
Valid filters are:\n\
    Tukey,a,q   - Tukey filter with parameters a & q\n\
    Hamming,a,b - Hamming filter with parameters a & b\n\
    Gauss,a     - Gaussian with width a (in k-space)\n"
};

int main(int argc, char **argv) {
    Eigen::initParallel();

    //******************************************************************************
    // Arguments / Usage
    //******************************************************************************
    bool verbose = false, debug = false;
    int complex_in = false, complex_out = false;
    shared_ptr<QI::FilterKernel> kernel = make_shared<QI::TukeyKernel>();
    string out_name;
    const struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, 0, 'v'},
        {"out", required_argument, 0, 'o'},
        {"filter", required_argument, 0, 'f'},
        {"debug", no_argument, 0, 'd'},
        {"complex_in", no_argument, &complex_in, true},
        {"complex_out", no_argument, &complex_out, true},
        {"threads", required_argument, 0, 'T'},
        {0, 0, 0, 0}
    };
    const char *short_options = "hvo:m:e:dT:";

    int indexptr = 0, c;
    while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
        switch (c) {
            case 'v': verbose = true; break;
            case 'o':
                out_name = optarg;
                if (verbose) cout << "Output filename will be: " << out_name << endl;
                break;
            case 'f': {
                stringstream f(optarg);
                kernel = QI::ReadKernel(f);
                if (verbose) cout << "Kernel is: " << *kernel << endl;
            } break;
            case 'd': debug = true; break;
            case 'T': itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg)); break;
            case 'h':
                cout << QI::GetVersion() << endl << usage << endl;
                return EXIT_SUCCESS;
            case 0: break; // Just a flag
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
    
    QI::SeriesXF::Pointer vols;
    if (complex_in) {
        if (verbose) cout << "Reading complex file: " << in_name << endl;
        vols = QI::ReadImage<QI::SeriesXF>(in_name);
    } else {
        if (verbose) cout << "Reading real file: " << in_name << endl;
        QI::SeriesF::Pointer rvols = QI::ReadImage<QI::SeriesF>(in_name);
        auto cast = itk::CastImageFilter<QI::SeriesF, QI::SeriesXF>::New();
        cast->SetInput(rvols);
        cast->Update();
        vols = cast->GetOutput();
        vols->DisconnectPipeline();
    }
    vols = run_pipeline(vols, verbose, debug, kernel);     
    if (complex_out) {
        if (verbose) cout << "Data is complex." << endl;
        QI::WriteImage(vols, out_name + QI::OutExt());
    } else {
        if (verbose) cout << "Data is real." << endl;
        QI::WriteMagnitudeImage(vols, out_name + QI::OutExt());
    }
    return EXIT_SUCCESS;
}
