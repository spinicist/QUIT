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

#include "Util.h"
#include "ImageIO.h"
#include "ApplyTypes.h"
#include "Args.h"

#include "itkRegionOfInterestImageFilter.h"
#include "itkStatisticsImageFilter.h"

//******************************************************************************
// Algorithm Subclasses
//******************************************************************************
typedef itk::ApplyAlgorithmFilter<QI::VectorVolumeXF, QI::VectorVolumeXF, QI::VectorVolumeF, QI::VolumeF> TApplyCombine;
class ComplexCombine : public TApplyCombine::Algorithm {
protected:
    const int m_insize = 0, m_coils = 0, m_outsize = 0;
    TOutput m_zero;
    std::vector<TConst> m_channel_phase;
public:
    ComplexCombine(const int insize, const int coils) :
        m_insize(insize), m_coils(coils), m_outsize(insize / coils)
    {
        TConst temp(m_coils);
        m_channel_phase.clear();
        m_channel_phase.push_back(temp);
        m_zero = TOutput(m_outsize);
        m_zero.Fill(std::complex<float>(0, 0));
        // QI_DB(m_insize)
        // QI_DB(m_coils)
        // QI_DB(m_outsize)
    }
    void setChannelPhases(const TConst &ph) {
        if (ph.Size() != m_coils) {
            QI_FAIL("Number of channel reference phases " << ph.Size() << " not equal to number of input coils " << m_coils);
        }
        m_channel_phase.clear();
        m_channel_phase.push_back(ph);
    }
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 1; }
    size_t numOutputs() const override { return 1; }
    size_t dataSize() const override { return m_insize; }
    size_t outputSize() const override { return m_outsize; }
    std::vector<TConst> defaultConsts() const override {
        return m_channel_phase;
    }
    TOutput zero() const override {
        return m_zero;
    }

    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TOutput &residual,
               TInput &resids, TIterations &its) const override
    {
        Eigen::Map<const Eigen::ArrayXXcf> in_data(inputs[0].GetDataPointer(), m_outsize, m_coils);
        Eigen::ArrayXXcd data = in_data.cast<std::complex<double>>();

        Eigen::Map<const Eigen::ArrayXf> in_phase(consts[0].GetDataPointer(), m_coils);
        Eigen::ArrayXd ph   = in_phase.cast<double>();
        Eigen::ArrayXcd correction(m_coils);
        correction.real() = ph.cos();
        correction.imag() = ph.sin();

        // Remove phase
        Eigen::ArrayXXcd phase_corrected = data.rowwise() / correction.transpose();
        // Then average
        Eigen::ArrayXcd averaged = phase_corrected.rowwise().sum() / m_coils;

        // QI_DBMAT(data.abs())
        // QI_DBMAT(data.arg())
        // QI_DBVEC(ph)
        // QI_DBVEC(correction.abs())
        // QI_DBVEC(correction.arg())
        // QI_DBMAT(phase_corrected.abs())
        // QI_DBMAT(phase_corrected.arg())
        // QI_DBVEC(averaged.abs())
        // QI_DBVEC(averaged.arg())

        for (auto i = 0; i < m_outsize; i++) {
            outputs[0][i] = averaged[i];
        }
        Eigen::ArrayXcf pcf = (Eigen::Map<Eigen::ArrayXcd>(phase_corrected.data(), m_insize)).cast<std::complex<float>>();
        for (auto i = 0; i < resids.Size(); i++) { // resids will be zero-length if not saving residuals
            resids[i] = pcf[i];
        }
        its = 1;
        return true;
    }
};

class ComplexVectorMeanFilter : public itk::ImageToImageFilter<QI::VectorVolumeXF, QI::VectorVolumeXF> {
public:
    /** Standard class typedefs. */
    typedef QI::VectorVolumeXF     TImage;

    typedef ComplexVectorMeanFilter            Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef itk::SmartPointer<Self>            Pointer;
    typedef typename TImage::RegionType        RegionType;
    typedef typename TImage::PixelType         PixelType;
    
    itkNewMacro(Self);
    itkTypeMacro(Self, Superclass);

    const PixelType &GetResult() const { return m_mean; } 
    void GenerateOutputInformation() ITK_OVERRIDE {
        Superclass::GenerateOutputInformation();
        auto op = this->GetOutput();
        op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
        op->Allocate();
    }

protected:
    PixelType m_mean;

    ComplexVectorMeanFilter() {
        this->SetNumberOfRequiredInputs(1);
    }
    ~ComplexVectorMeanFilter() {}

    void GenerateData() ITK_OVERRIDE {
        typename TImage::ConstPointer input = this->GetInput();
        auto region = input->GetLargestPossibleRegion();
        auto N = region.GetNumberOfPixels();
        m_mean = PixelType(input->GetNumberOfComponentsPerPixel());
        m_mean.Fill(std::complex<float>(0.,0.));
        itk::ImageRegionConstIterator<TImage> imageIt(input,region);
        imageIt.GoToBegin();
        while(!imageIt.IsAtEnd()) {
            m_mean += imageIt.Get();
            ++imageIt;
        }
        m_mean /= N;
    }

private:
    ComplexVectorMeanFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

int main(int argc, char **argv) {
    args::ArgumentParser parser(
        "Combine multiple coil images into a single image.\n"
        "Default method is that of Hammond, 10.1016/j.neuroimage.2007.10.037\n"
        "If the COMPOSER option is specified, see 10.1002/mrm.26093\n"
        "http://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "INPUT_FILE", "Input file to coil-combine");
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> region_arg(parser, "REGION", "Region to average phase for Hammond method, default is 8x8x8 cube at center", {'r', "region"});
    args::ValueFlag<std::string> ser_path(parser, "COMPOSER", "Short Echo Time reference file for COMPOSER method", {'c', "composer"});
    args::ValueFlag<int> coils_arg(parser, "COILS", "Number of coils (default is number of volumes)", {'C', "coils"});
    args::Flag     save_corrected(parser, "SAVE COILS", "Save the individual coil images after phase correction", {'s', "save"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv, verbose);

    if (verbose) std::cout << "Reading input image: " << QI::CheckPos(input_path) << std::endl;
    auto input_image = QI::ReadVectorImage<std::complex<float>>(QI::CheckPos(input_path));
    const auto sz = input_image->GetNumberOfComponentsPerPixel();
    const auto ncoils = coils_arg ? coils_arg.Get() : sz;
    auto combine = std::make_shared<ComplexCombine>(sz, ncoils);
    auto apply = TApplyCombine::New();
    apply->SetAlgorithm(combine);
    apply->SetInput(0, input_image);
    apply->SetOutputAllResiduals(save_corrected);
    apply->SetVerbose(verbose);
    apply->SetPoolsize(threads.Get());
    apply->SetSplitsPerThread(threads.Get());
    if (subregion) apply->SetSubregion(QI::RegionArg(subregion.Get()));
    if (ser_path) {
        if (verbose) std::cout << "Reading COMPOSER reference image: " << ser_path.Get() << std::endl;
        auto ser_image = QI::ReadVectorImage(ser_path.Get());
        if (ser_image->GetNumberOfComponentsPerPixel() != ncoils) {
            QI_FAIL("Number of coil reference images does not match number of coils in data");
        }
        apply->SetConst(0, ser_image);
    } else {
        // Fall back to Hammond Method
        if (verbose) std::cout << "Using Hammond method" << std::endl;
        QI::VectorVolumeXF::RegionType region;
        if (region_arg) {
            region = QI::RegionArg<QI::VolumeF::RegionType>(region_arg.Get());
        } else {
            auto size = input_image->GetLargestPossibleRegion().GetSize();
            itk::Index<3> index;
            for (auto i = 0; i < 3; i++) {
                index[i] = size[i] / 2 - 4;
            }
            region.GetModifiableIndex() = index;
            region.GetModifiableSize() = {{8, 8, 8}};
        }
        if (verbose) std::cout << "Reference region is:\n" << region << std::endl;
        auto roi = itk::RegionOfInterestImageFilter<QI::VectorVolumeXF, QI::VectorVolumeXF>::New();
        roi->SetRegionOfInterest(region);
        roi->SetInput(input_image);
        auto mean_filter = ComplexVectorMeanFilter::New();
        mean_filter->SetInput(roi->GetOutput());
        mean_filter->Update();
        auto roi_mean = mean_filter->GetResult();
        if (verbose) std::cout << "Mean values: " << roi_mean << std::endl;
        itk::VariableLengthVector<float> phase(sz);
        for (auto i = 0; i < sz; i++) {
            phase[i] = std::arg(roi_mean[i]);
        }
        if (verbose) std::cout << "Mean phase: " << phase << std::endl;
        combine->setChannelPhases(phase);
    }
    if (verbose) std::cout << "Correcting phase & combining" << std::endl;
    apply->Update();
    const std::string out_name = (outarg ? outarg.Get() : QI::StripExt(input_path.Get())) + "_combined" + QI::OutExt();
    if (verbose) std::cout << "Writing output file " << out_name << std::endl;
    QI::WriteVectorImage(apply->GetOutput(0), out_name);
    if (save_corrected) {
        const std::string out_name = (outarg ? outarg.Get() : QI::StripExt(input_path.Get())) + "_corrected" + QI::OutExt();
        if (verbose) std::cout << "Writing corrected coil file " << out_name << std::endl;
        QI::WriteVectorImage(apply->GetAllResidualsOutput(), out_name);
    }
    return EXIT_SUCCESS;
}

