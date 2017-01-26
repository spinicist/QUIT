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

template<typename TImage>
class KernelSource : public ImageSource<TImage> {
public:
    typedef TImage                 ImageType;
    typedef KernelSource           Self;
    typedef ImageSource<ImageType> Superclass;
    typedef SmartPointer<Self>     Pointer;
    typedef typename ImageType::RegionType     RegionType;
    typedef typename ImageType::IndexType      IndexType;
    typedef typename ImageType::SizeType       SizeType;
    typedef typename ImageType::SpacingType    SpacingType;
    typedef typename ImageType::DirectionType  DirectionType;
    typedef typename ImageType::PointType      PointType;

protected:
    KernelSource(){}
    ~KernelSource(){}
    RegionType    m_Region;
    SpacingType   m_Spacing;
    DirectionType m_Direction;
    PointType     m_Origin;
    std::vector<std::shared_ptr<QI::FilterKernel>> m_kernels;

public:
    itkNewMacro(Self);
    itkTypeMacro(Self, ImageSource);
    void SetKernels(const std::vector<std::shared_ptr<QI::FilterKernel>> &k) { m_kernels = k; }
    itkSetMacro(Region, RegionType);
    itkSetMacro(Spacing, SpacingType);
    itkSetMacro(Direction, DirectionType);
    itkSetMacro(Origin, PointType);

protected:
    virtual void GenerateOutputInformation() ITK_OVERRIDE {
        auto output = this->GetOutput(0);
        output->SetLargestPossibleRegion(m_Region);
        output->SetSpacing(m_Spacing);
        output->SetDirection(m_Direction);
        output->SetOrigin(m_Origin);
    }

    virtual void ThreadedGenerateData(const RegionType &region, ThreadIdType threadId) ITK_OVERRIDE {
        const auto startIndex = m_Region.GetIndex();
        const Eigen::Array3d sz{m_Region.GetSize()[0], m_Region.GetSize()[1], m_Region.GetSize()[2]};
        const Eigen::Array3d hsz = sz / 2;
        const Eigen::Array3d sp{m_Spacing[0], m_Spacing[1], m_Spacing[2]};
        itk::ImageRegionIterator<ImageType> outIt(this->GetOutput(), m_Region);
        outIt.GoToBegin();
        while(!outIt.IsAtEnd()) {
            const auto I = outIt.GetIndex() - startIndex; // Might be padded to a negative start
            const double x = fmod(static_cast<double>(I[0]) + hsz[0], sz[0]) - hsz[0];
            const double y = fmod(static_cast<double>(I[1]) + hsz[1], sz[1]) - hsz[1];
            const double z = fmod(static_cast<double>(I[2]) + hsz[2], sz[2]) - hsz[2];
            const Eigen::Array3d p{x, y, z};
            typename ImageType::PixelType val = 1;
            for (const auto &kernel : m_kernels) {
                val *= kernel->value(p, hsz, sp);
            }
            outIt.Set(val);
            ++outIt;
        }
    }

private:
    KernelSource(const Self &);
    void operator=(const Self &);
};

} // End namespace itk

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
    --save_kernel, -k    : Save all pipeline steps.\n\
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
    bool verbose = false, save_kernel = false;
    int complex_in = false, complex_out = false;
    std::vector<shared_ptr<QI::FilterKernel>> kernels;
    string out_name;
    const struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, 0, 'v'},
        {"out", required_argument, 0, 'o'},
        {"filter", required_argument, 0, 'f'},
        {"save_kernel", no_argument, 0, 'k'},
        {"complex_in", no_argument, &complex_in, true},
        {"complex_out", no_argument, &complex_out, true},
        {"threads", required_argument, 0, 'T'},
        {0, 0, 0, 0}
    };
    const char *short_options = "hvo:m:e:kT:";

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
                kernels.push_back(QI::ReadKernel(f));
                if (verbose) cout << "Read kernel: " << *kernels.back() << endl;
            } break;
            case 'k': save_kernel = true; break;
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

    /*
     * Main calculation starts hear.
     * First a lot of typedefs.
     */
    typedef itk::ExtractImageFilter<QI::SeriesXF, QI::VolumeXD> TExtract;
    typedef itk::TileImageFilter<QI::VolumeXF, QI::SeriesXF>    TTile;
    typedef itk::CastImageFilter<QI::VolumeXD, QI::VolumeXF>    TCast;
    typedef itk::FFTPadImageFilter<QI::VolumeXD>                TPad;
    typedef itk::ComplexToComplexFFTImageFilter<QI::VolumeXD>   TFFT;
    typedef itk::MultiplyImageFilter<QI::VolumeXD, QI::VolumeD, QI::VolumeXD> TMult;
    typedef itk::KernelSource<QI::VolumeD>                      TKernel;
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
    auto tkernel = TKernel::New();
    auto mult    = TMult::New();

    // inverse is declared in the loop due to a weird bug
    auto unpadder = TUnpad::New();
    auto tile   = TTile::New();
    
    extract->SetInput(vols);
    extract->SetDirectionCollapseToSubmatrix();
    pad->SetSizeGreatestPrimeFactor(5); // This is the largest the VNL FFT supports
    pad->SetInput(extract->GetOutput());
    forward->SetInput(pad->GetOutput());
    if (kernels.size() == 0) {
        kernels.push_back(std::make_shared<QI::TukeyKernel>());
    }
    tkernel->SetKernels(kernels);
    mult->SetInput1(forward->GetOutput());
    mult->SetInput2(tkernel->GetOutput());
    unpadder->SetDirectionCollapseToSubmatrix();
    unpadder->SetExtractionRegion(unpad_region);
    tile->SetLayout(layout);
    
    for (int i = 0; i < nvols; i++) {
        region.GetModifiableIndex()[3] = i;
        if (verbose) cout << "Processing volume " << i << endl;
        
        extract->SetExtractionRegion(region);
        extract->Update();
        pad->Update();
        if (i == 0) { // For first image we need to update the kernel
            tkernel->SetRegion(pad->GetOutput()->GetLargestPossibleRegion());
            tkernel->SetSpacing(pad->GetOutput()->GetSpacing());
            tkernel->SetOrigin(pad->GetOutput()->GetOrigin());
            tkernel->SetDirection(pad->GetOutput()->GetDirection());
            tkernel->Update();
            if (verbose) std::cout << "Built kernel" << std::endl;
        }
        mult->Update();
        auto inverse = TFFT::New();
        inverse->SetTransformDirection(TFFT::INVERSE);
        inverse->SetInput(mult->GetOutput());
        unpadder->SetInput(inverse->GetOutput());
        inverse->Update(); // Need a separate update to avoid region bug
        unpadder->Update();
        QI::VolumeXF::Pointer v = unpadder->GetOutput();
        tile->SetInput(i, v);
        v->DisconnectPipeline();
    }
    if (verbose) cout << "Finished." << endl;
    tile->Update();
    auto dir = vols->GetDirection();
    auto spc = vols->GetSpacing();
    vols = tile->GetOutput();
    vols->SetDirection(dir);
    vols->SetSpacing(spc);
    vols->DisconnectPipeline();

    if (out_name == "")
        out_name = in_name.substr(0, in_name.find("."));
    const std::string out_path = out_name + "_filtered" + QI::OutExt();
    if (complex_out) {
        if (verbose) cout << "Saving complex output file: " << out_path << endl;
        QI::WriteImage(vols, out_path);
    } else {
        if (verbose) cout << "Saving real output file: " << out_path << endl;
        QI::WriteMagnitudeImage(vols, out_path);
    }
    if (save_kernel) {
        const string kernel_path = out_name + "_kernel" + QI::OutExt();
        if (verbose) cout << "Saving filter kernel to: " << kernel_path << endl;
        QI::WriteImage(tkernel->GetOutput(), kernel_path);
    }
    return EXIT_SUCCESS;
}
