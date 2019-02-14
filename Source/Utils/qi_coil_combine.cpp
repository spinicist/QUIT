/*
 *  qicomposer.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <array>

#include <Eigen/Core>

#include "Args.h"
#include "ImageIO.h"
#include "Util.h"

#include "itkRegionOfInterestImageFilter.h"
#include "itkStatisticsImageFilter.h"

class CoilCombineFilter : public itk::ImageToImageFilter<QI::VectorVolumeXF, QI::VectorVolumeXF> {
  public:
    /** Standard class typedefs. */
    typedef QI::VectorVolumeXF TImage;

    typedef CoilCombineFilter                  Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef itk::SmartPointer<Self>            Pointer;
    typedef typename TImage::RegionType        RegionType;
    typedef typename TImage::PixelType         PixelType;

    itkNewMacro(Self);
    itkTypeMacro(Self, Superclass);

    void GenerateOutputInformation() ITK_OVERRIDE {
        Superclass::GenerateOutputInformation();
        auto input        = this->GetInput(0);
        auto region       = input->GetLargestPossibleRegion();
        auto spacing      = input->GetSpacing();
        auto origin       = input->GetOrigin();
        auto direction    = input->GetDirection();
        m_coils           = this->GetInput(1)->GetNumberOfComponentsPerPixel();
        m_images_per_coil = input->GetNumberOfComponentsPerPixel() / m_coils;
        auto op           = this->GetOutput(0);
        op->SetRegions(region);
        op->SetSpacing(spacing);
        op->SetOrigin(origin);
        op->SetDirection(direction);
        op->SetNumberOfComponentsPerPixel(m_images_per_coil);
        op->Allocate(true);
    }

  protected:
    int m_coils = 1, m_images_per_coil = 1;

    CoilCombineFilter() { this->SetNumberOfRequiredInputs(2); }
    ~CoilCombineFilter() {}

    void DynamicThreadedGenerateData(const RegionType &region) ITK_OVERRIDE {
        const auto                                        input_image  = this->GetInput(0);
        const auto                                        ref_image    = this->GetInput(1);
        auto                                              output_image = this->GetOutput();
        itk::ImageRegionConstIterator<QI::VectorVolumeXF> input_iter(input_image, region);
        itk::ImageRegionConstIterator<QI::VectorVolumeXF> ref_iter(ref_image, region);
        itk::ImageRegionIterator<QI::VectorVolumeXF>      output_iter(output_image, region);
        input_iter.GoToBegin();
        ref_iter.GoToBegin();
        output_iter.GoToBegin();

        while (!input_iter.IsAtEnd()) {
            const auto input_vec = input_iter.Get();
            const auto ref_vec   = ref_iter.Get();

            Eigen::Map<const Eigen::ArrayXXcf> data(input_vec.GetDataPointer(), m_images_per_coil,
                                                    m_coils);
            Eigen::Map<const Eigen::ArrayXcf>  ref(ref_vec.GetDataPointer(), m_coils);
            const Eigen::ArrayXcf              correction = ref / ref.abs();
            // Remove phase
            const Eigen::ArrayXXcf phase_corrected = data.rowwise() / correction.transpose();
            // Then average
            const Eigen::ArrayXcf averaged = phase_corrected.rowwise().sum() / m_coils;

            for (int i = 0; i < m_images_per_coil; i++) {
                output_iter.Get()[i] = averaged[i];
            }

            ++input_iter;
            ++ref_iter;
            ++output_iter;
        }
    }

  private:
    CoilCombineFilter(const Self &); // purposely not implemented
    void operator=(const Self &);    // purposely not implemented
};

class HammondCombineFilter
    : public itk::ImageToImageFilter<QI::VectorVolumeXF, QI::VectorVolumeXF> {
  public:
    /** Standard class typedefs. */
    typedef QI::VectorVolumeXF TImage;

    typedef HammondCombineFilter               Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef itk::SmartPointer<Self>            Pointer;
    typedef typename TImage::RegionType        RegionType;
    typedef typename TImage::PixelType         PixelType;

    itkNewMacro(Self);
    itkTypeMacro(Self, Superclass);

    void GenerateOutputInformation() ITK_OVERRIDE {
        Superclass::GenerateOutputInformation();
        auto input     = this->GetInput(0);
        auto region    = input->GetLargestPossibleRegion();
        auto spacing   = input->GetSpacing();
        auto origin    = input->GetOrigin();
        auto direction = input->GetDirection();

        m_coils           = m_hammond_ref.rows();
        m_images_per_coil = input->GetNumberOfComponentsPerPixel() / m_coils;
        auto op           = this->GetOutput(0);
        op->SetRegions(region);
        op->SetSpacing(spacing);
        op->SetOrigin(origin);
        op->SetDirection(direction);
        op->SetNumberOfComponentsPerPixel(m_images_per_coil);
        op->Allocate(true);
    }

    void SetHammondRef(const Eigen::ArrayXcf &h) { m_hammond_ref = h / h.abs(); }

  protected:
    int             m_coils = 1, m_images_per_coil = 1;
    Eigen::ArrayXcf m_hammond_ref;

    HammondCombineFilter() { this->SetNumberOfRequiredInputs(1); }
    ~HammondCombineFilter() {}

    void DynamicThreadedGenerateData(const RegionType &region) ITK_OVERRIDE {
        const auto                                        input_image  = this->GetInput(0);
        auto                                              output_image = this->GetOutput();
        itk::ImageRegionConstIterator<QI::VectorVolumeXF> input_iter(input_image, region);
        itk::ImageRegionIterator<QI::VectorVolumeXF>      output_iter(output_image, region);
        input_iter.GoToBegin();
        output_iter.GoToBegin();

        while (!input_iter.IsAtEnd()) {
            const auto input_vec = input_iter.Get();

            Eigen::Map<const Eigen::ArrayXXcf> data(input_vec.GetDataPointer(), m_images_per_coil,
                                                    m_coils);
            // Remove phase
            const Eigen::ArrayXXcf phase_corrected = data.rowwise() / m_hammond_ref.transpose();
            // Then average
            const Eigen::ArrayXcf averaged = phase_corrected.rowwise().sum() / m_coils;

            for (int i = 0; i < m_images_per_coil; i++) {
                output_iter.Get()[i] = averaged[i];
            }

            ++input_iter;
            ++output_iter;
        }
    }

  private:
    HammondCombineFilter(const Self &); // purposely not implemented
    void operator=(const Self &);       // purposely not implemented
};

class ComplexVectorMeanFilter
    : public itk::ImageToImageFilter<QI::VectorVolumeXF, QI::VectorVolumeXF> {
  public:
    /** Standard class typedefs. */
    typedef QI::VectorVolumeXF TImage;

    typedef ComplexVectorMeanFilter            Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef itk::SmartPointer<Self>            Pointer;
    typedef typename TImage::RegionType        RegionType;
    typedef typename TImage::PixelType         PixelType;

    itkNewMacro(Self);
    itkTypeMacro(Self, Superclass);

    const PixelType &GetResult() const { return m_mean; }
    void             GenerateOutputInformation() ITK_OVERRIDE {
        Superclass::GenerateOutputInformation();
        auto op = this->GetOutput();
        op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
        op->Allocate();
    }

  protected:
    PixelType m_mean;

    ComplexVectorMeanFilter() { this->SetNumberOfRequiredInputs(1); }
    ~ComplexVectorMeanFilter() {}

    void GenerateData() ITK_OVERRIDE {
        typename TImage::ConstPointer input  = this->GetInput();
        auto                          region = input->GetLargestPossibleRegion();
        auto                          N      = region.GetNumberOfPixels();
        m_mean                               = PixelType(input->GetNumberOfComponentsPerPixel());
        m_mean.Fill(std::complex<float>(0., 0.));
        itk::ImageRegionConstIterator<TImage> imageIt(input, region);
        imageIt.GoToBegin();
        while (!imageIt.IsAtEnd()) {
            m_mean += imageIt.Get();
            ++imageIt;
        }
        m_mean /= N;
    }

  private:
    ComplexVectorMeanFilter(const Self &); // purposely not implemented
    void operator=(const Self &);          // purposely not implemented
};

int main(int argc, char **argv) {
    args::ArgumentParser parser(
        "Combine multiple coil images into a single image.\n"
        "Default method is that of Hammond, 10.1016/j.neuroimage.2007.10.037\n"
        "If the COMPOSER option is specified, see 10.1002/mrm.26093\n"
        "http://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "INPUT_FILE", "Input file to coil-combine");
    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames",
                                        {'o', "out"});
    args::ValueFlag<std::string> region_arg(
        parser, "REGION",
        "Region to average phase for Hammond method, default is 8x8x8 cube at center",
        {'r', "region"});
    args::ValueFlag<std::string> ser_path(parser, "COMPOSER",
                                          "Short Echo Time reference file for COMPOSER method",
                                          {'c', "composer"});
    args::ValueFlag<int>         coils_arg(
        parser, "COILS", "Number of coils for Hammond method (default is number of volumes)",
        {'C', "coils"}, -1);
    args::ValueFlag<int> ref_vol(parser, "VOLUME",
                                 "Volume to use as reference for Hammond method (default is 1)",
                                 {'V', "vol"}, 1);
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    auto input_image = QI::ReadImage<QI::VectorVolumeXF>(QI::CheckPos(input_path), verbose);

    QI::VectorVolumeXF::Pointer output = ITK_NULLPTR;
    if (ser_path) {
        QI::Log(verbose, "Reading COMPOSER reference image: {}", ser_path.Get());
        auto ser_image = QI::ReadImage<QI::VectorVolumeXF>(ser_path.Get(), verbose);
        auto combine   = CoilCombineFilter::New();
        combine->SetInput(input_image);
        combine->SetInput(1, ser_image);
        QI::Log(verbose, "Applying COMPOSER");
        combine->Update();
        output = combine->GetOutput();
    } else {
        // Fall back to Hammond Method
        QI::Log(verbose, "Using Hammond method");
        QI::VectorVolumeXF::RegionType region;
        if (region_arg) {
            region = QI::RegionFromString<QI::VolumeF::RegionType>(region_arg.Get());
        } else {
            auto          size = input_image->GetLargestPossibleRegion().GetSize();
            itk::Index<3> index;
            for (auto i = 0; i < 3; i++) {
                index[i] = size[i] / 2 - 4;
            }
            region.GetModifiableIndex() = index;
            region.GetModifiableSize()  = {{8, 8, 8}};
        }
        QI::Log(verbose, "Reference region is:\n{}", region);
        auto roi = itk::RegionOfInterestImageFilter<QI::VectorVolumeXF, QI::VectorVolumeXF>::New();
        roi->SetRegionOfInterest(region);
        roi->SetInput(input_image);
        auto mean_filter = ComplexVectorMeanFilter::New();
        mean_filter->SetInput(roi->GetOutput());
        mean_filter->Update();
        auto roi_mean = mean_filter->GetResult();
        QI::Log(verbose, "Mean values: {}", roi_mean);

        Eigen::Map<const Eigen::ArrayXcf> mean(roi_mean.GetDataPointer(), roi_mean.Size(), 1);
        Eigen::ArrayXcf                   hammond_ref(coils_arg.Get());
        const int                         images_per_coil = mean.rows() / coils_arg.Get();
        for (int i = 0; i < coils_arg.Get(); i++) {
            hammond_ref(i) = mean(ref_vol.Get() + i * images_per_coil);
        }
        QI::Log(verbose, "Hammond ref: {}", hammond_ref.transpose());
        auto hammond = HammondCombineFilter::New();
        hammond->SetInput(input_image);
        hammond->SetHammondRef(hammond_ref);
        QI::Log(verbose, "Applying Hammond method");
        hammond->Update();
        output = hammond->GetOutput();
        // fit_filter->SetFixed(0, phase);
    }
    QI::Log(verbose, "Correcting phase & combining");
    const std::string out_name =
        (outarg ? outarg.Get() : QI::StripExt(input_path.Get())) + "_combined" + QI::OutExt();
    QI::WriteImage(output, out_name, verbose);
    return EXIT_SUCCESS;
}
