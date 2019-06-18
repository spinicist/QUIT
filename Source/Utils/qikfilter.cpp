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

#include "Eigen/Core"
#include <sstream>

#include "itkCastImageFilter.h"
#include "itkComplexToComplexFFTImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkConstantPadImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkFFTPadImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageSource.h"
#include "itkInverseFFTImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkTileImageFilter.h"

#include "Args.h"
#include "ImageIO.h"
#include "ImageTypes.h"
#include "Kernels.h"
#include "Util.h"

using namespace Eigen;

namespace itk {

template <typename TImage> class KernelSource : public ImageSource<TImage> {
  public:
    typedef TImage                            ImageType;
    typedef KernelSource                      Self;
    typedef ImageSource<ImageType>            Superclass;
    typedef SmartPointer<Self>                Pointer;
    typedef typename ImageType::RegionType    RegionType;
    typedef typename ImageType::IndexType     IndexType;
    typedef typename ImageType::SizeType      SizeType;
    typedef typename ImageType::SpacingType   SpacingType;
    typedef typename ImageType::DirectionType DirectionType;
    typedef typename ImageType::PointType     PointType;

  protected:
    KernelSource() { this->DynamicMultiThreadingOff(); }
    ~KernelSource() {}
    RegionType                                     m_Region;
    SpacingType                                    m_Spacing;
    DirectionType                                  m_Direction;
    PointType                                      m_Origin;
    std::vector<std::shared_ptr<QI::FilterKernel>> m_kernels;
    bool                                           m_Highpass = false;

  public:
    itkNewMacro(Self);
    itkTypeMacro(Self, ImageSource);
    void SetKernel(const std::shared_ptr<QI::FilterKernel> &k) {
        m_kernels.clear();
        m_kernels.push_back(k);
        this->Modified();
    }
    void SetKernels(const std::vector<std::shared_ptr<QI::FilterKernel>> &k) {
        m_kernels = k;
        this->Modified();
    }
    itkSetMacro(Region, RegionType);
    itkSetMacro(Spacing, SpacingType);
    itkSetMacro(Direction, DirectionType);
    itkSetMacro(Origin, PointType);
    itkSetMacro(Highpass, bool);

  protected:
    void GenerateOutputInformation() ITK_OVERRIDE {
        auto output = this->GetOutput(0);
        output->SetLargestPossibleRegion(m_Region);
        output->SetSpacing(m_Spacing);
        output->SetDirection(m_Direction);
        output->SetOrigin(m_Origin);
    }

    void DynamicThreadedGenerateData(const RegionType &outputRegion) ITK_OVERRIDE {
        const auto                          startIndex = m_Region.GetIndex();
        const Eigen::Array3d                sz{static_cast<double>(m_Region.GetSize()[0]),
                                static_cast<double>(m_Region.GetSize()[1]),
                                static_cast<double>(m_Region.GetSize()[2])};
        const Eigen::Array3d                hsz = sz / 2;
        const Eigen::Array3d                sp{m_Spacing[0], m_Spacing[1], m_Spacing[2]};
        itk::ImageRegionIterator<ImageType> outIt(this->GetOutput(), outputRegion);
        outIt.GoToBegin();
        while (!outIt.IsAtEnd()) {
            const auto   I = outIt.GetIndex() - startIndex; // Might be padded to a negative start
            const double x = fmod(static_cast<double>(I[0]) + hsz[0], sz[0]) - hsz[0];
            const double y = fmod(static_cast<double>(I[1]) + hsz[1], sz[1]) - hsz[1];
            const double z = fmod(static_cast<double>(I[2]) + hsz[2], sz[2]) - hsz[2];
            const Eigen::Array3d          p{x, y, z};
            typename ImageType::PixelType val = 1;
            for (const auto &kernel : m_kernels) {
                val *= kernel->value(p, hsz, sp);
            }
            if (m_Highpass) {
                val = 1 - val;
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

int main(int argc, char **argv) {
    Eigen::initParallel();

    args::ArgumentParser parser("smooths images in k-space\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> in_path(parser, "INPUT", "Input file.");
    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string> out_prefix(
        parser, "OUTPREFIX", "Change output prefix (default input filename)", {'o', "out"});
    args::ValueFlag<int> zero_padding(
        parser, "ZEROPAD", "Zero-pad volume by N voxels in each direction", {'z', "zero_pad"}, 0);
    args::Flag complex_in(parser, "COMPLEX_IN", "Input data is complex", {"complex_in"});
    args::Flag complex_out(parser, "COMPLEX_OUT", "Write complex output", {"complex_out"});
    args::Flag save_kernel(parser, "KERNEL", "Save kernels as images", {"save_kernel"});
    args::Flag save_kspace(
        parser, "KSPACE", "Save k-space before & after filtering", {"save_kspace"});
    args::Flag highpass(parser, "HIGHPASS", "Use a high-pass, not a low-pass filter", {"highpass"});
    args::Flag filter_per_volume(parser,
                                 "FILTER_PER_VOL",
                                 "Instead of concatenating multiple filters, use one per volume",
                                 {"filter_per_volume"});
    args::ValueFlagList<std::string> filters(
        parser, "FILTER", "Specify a filter to use (can be multiple)", {'f', "filter"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    std::vector<std::shared_ptr<QI::FilterKernel>> kernels;
    if (filters) {
        for (const auto &f : filters.Get()) {
            kernels.push_back(QI::ReadKernel(f));
            QI::Log(verbose, "Read kernel: {}", *(kernels.back()));
        }
    } else {
        kernels.push_back(std::make_shared<QI::TukeyKernel>());
    }

    QI::SeriesXF::Pointer vols;
    if (complex_in) {
        QI::Log(verbose, "Reading complex file: {}", QI::CheckPos(in_path));
        vols = QI::ReadImage<QI::SeriesXF>(QI::CheckPos(in_path), verbose);
    } else {
        QI::Log(verbose, "Reading real file: {}", QI::CheckPos(in_path));
        QI::SeriesF::Pointer rvols = QI::ReadImage<QI::SeriesF>(QI::CheckPos(in_path), verbose);
        auto                 cast  = itk::CastImageFilter<QI::SeriesF, QI::SeriesXF>::New();
        cast->SetInput(rvols);
        cast->Update();
        vols = cast->GetOutput();
        vols->DisconnectPipeline();
    }
    const std::string out_base = out_prefix ? out_prefix.Get() : QI::Basename(in_path.Get());
    /*
     * Main calculation starts hear.
     * First a lot of typedefs.
     */
    typedef itk::ExtractImageFilter<QI::SeriesXF, QI::VolumeXD>               TExtract;
    typedef itk::TileImageFilter<QI::VolumeXF, QI::SeriesXF>                  TTile;
    typedef itk::ConstantPadImageFilter<QI::VolumeXD, QI::VolumeXD>           TZeroPad;
    typedef itk::FFTPadImageFilter<QI::VolumeXD>                              TFFTPad;
    typedef itk::FFTShiftImageFilter<QI::VolumeXD, QI::VolumeXD>              TFFTShift;
    typedef itk::ComplexToComplexFFTImageFilter<QI::VolumeXD>                 TFFT;
    typedef itk::MultiplyImageFilter<QI::VolumeXD, QI::VolumeD, QI::VolumeXD> TMult;
    typedef itk::KernelSource<QI::VolumeD>                                    TKernel;
    typedef itk::ExtractImageFilter<QI::VolumeXD, QI::VolumeXF>               TUnpad;

    auto         region = vols->GetLargestPossibleRegion();
    const size_t nvols  = region.GetSize()[3]; // Save for the loop
    if (filter_per_volume && nvols != kernels.size()) {
        QI::Fail(
            "Number of volumes ({}) and kernels ({}) do not match for filter_per_volume option",
            nvols,
            kernels.size());
    }
    region.GetModifiableSize()[3] = 0;
    QI::VolumeXD::RegionType unpad_region;
    for (int i = 0; i < 3; i++) {
        unpad_region.GetModifiableSize()[i] = region.GetSize()[i];
    }

    auto tkernel = TKernel::New();
    if (!filter_per_volume) {
        QI::Info(verbose, "Kernels:");
        for (auto const &k : kernels) {
            QI::Info(verbose, "{}", *k);
        }
        tkernel->SetKernels(kernels);
    }
    if (highpass) {
        tkernel->SetHighpass(highpass);
        QI::Log(verbose, "Set highpass filter");
    }

    auto                             tile = TTile::New();
    itk::FixedArray<unsigned int, 4> layout;
    layout[0] = layout[1] = layout[2] = 1;
    layout[3]                         = nvols;
    tile->SetLayout(layout);

    for (size_t i = 0; i < nvols; i++) {
        region.GetModifiableIndex()[3] = i;
        QI::Log(verbose, "Processing volume {}", i);

        auto extract = TExtract::New();
        extract->SetInput(vols);
        extract->SetDirectionCollapseToSubmatrix();
        extract->SetExtractionRegion(region);
        auto zero_pad = TZeroPad::New();
        auto fft_pad  = TFFTPad::New();
        fft_pad->SetSizeGreatestPrimeFactor(5); // This is the largest the VNL FFT supports
        if (zero_padding > 0) {
            QI::VolumeXD::SizeType padding;
            padding.Fill(zero_padding);
            zero_pad->SetInput(extract->GetOutput());
            zero_pad->SetPadLowerBound(padding);
            zero_pad->SetPadUpperBound(padding);
            zero_pad->SetConstant(0);
            fft_pad->SetInput(zero_pad->GetOutput());
        } else {
            fft_pad->SetInput(extract->GetOutput());
        }
        fft_pad->Update(); // Need to know the size of this to set up the kernel properly
        QI::Log(verbose,
                "After FFT padding size is: {}",
                fft_pad->GetOutput()->GetLargestPossibleRegion().GetSize());
        if (i == 0) {
            tkernel->SetRegion(fft_pad->GetOutput()->GetLargestPossibleRegion());
            tkernel->SetSpacing(fft_pad->GetOutput()->GetSpacing());
            tkernel->SetOrigin(fft_pad->GetOutput()->GetOrigin());
            tkernel->SetDirection(fft_pad->GetOutput()->GetDirection());
            tkernel->Update();
            QI::Log(verbose,
                    "Created kernel filter, size is: {}",
                    tkernel->GetOutput()->GetLargestPossibleRegion().GetSize());
        }
        if (filter_per_volume) {
            tkernel->SetKernel(kernels.at(i));
            QI::Log(verbose, "Setting kernel to: {}", *kernels.at(i));
        }
        auto forward = TFFT::New();
        forward->SetInput(fft_pad->GetOutput());
        forward->Update();

        auto mult = TMult::New();
        mult->SetInput1(forward->GetOutput());
        mult->SetInput2(tkernel->GetOutput());
        mult->Update();

        auto inverse = TFFT::New();
        inverse->SetTransformDirection(TFFT::INVERSE);
        inverse->SetInput(mult->GetOutput());
        QI::Log(verbose, "Filtering k-space");
        inverse->Update(); // If we don't have this update, next step fails

        auto unpadder = TUnpad::New();
        unpadder->SetDirectionCollapseToSubmatrix();
        unpadder->SetExtractionRegion(unpad_region);
        unpadder->SetInput(inverse->GetOutput());
        QI::Log(verbose, "Unpadding & tiling output");
        unpadder->Update();

        QI::VolumeXF::Pointer v = unpadder->GetOutput();
        tile->SetInput(i, v);
        v->DisconnectPipeline();

        if (save_kspace) {
            auto shift_filter = TFFTShift::New();
            auto cast_filter  = itk::CastImageFilter<QI::VolumeXD, QI::VolumeXF>::New();
            shift_filter->SetInput(forward->GetOutput());
            cast_filter->SetInput(shift_filter->GetOutput());
            cast_filter->Update();
            QI::WriteMagnitudeImage(
                cast_filter->GetOutput(), out_base + "_kspace_before" + QI::OutExt(), verbose);
            shift_filter->SetInput(mult->GetOutput());
            cast_filter->Update();
            QI::WriteMagnitudeImage(
                cast_filter->GetOutput(), out_base + "_kspace_after" + QI::OutExt(), verbose);
        }
    }
    QI::Log(verbose, "Finished.");
    tile->Update();
    auto dir = vols->GetDirection();
    auto spc = vols->GetSpacing();
    vols     = tile->GetOutput();
    vols->SetDirection(dir);
    vols->SetSpacing(spc);
    vols->DisconnectPipeline();

    const std::string out_path = out_base + "_filtered" + QI::OutExt();
    if (complex_out) {
        QI::WriteImage(vols, out_path, verbose);
    } else {
        QI::WriteMagnitudeImage(vols, out_path, verbose);
    }
    if (save_kernel) {
        auto shift_filter = itk::FFTShiftImageFilter<QI::VolumeD, QI::VolumeD>::New();
        shift_filter->SetInput(tkernel->GetOutput());
        auto cast_filter = itk::CastImageFilter<QI::VolumeD, QI::VolumeF>::New();
        cast_filter->SetInput(shift_filter->GetOutput());
        cast_filter->Update();
        QI::WriteImage(cast_filter->GetOutput(), out_base + "_kernel" + QI::OutExt(), verbose);
    }
    return EXIT_SUCCESS;
}
