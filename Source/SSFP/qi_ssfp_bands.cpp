/*
 *  qissfpbands.cpp
 *
 *  Created by Tobias Wood on 14/03/2014.
 *  Copyright (c) 2014 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>

#include "itkUnaryFunctorImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"

#include "Util.h"
#include "Args.h"
#include "Banding.h"
#include "ImageIO.h"

namespace itk {

class MinEnergyFilter : public ImageToImageFilter<QI::VectorVolumeXF, QI::VectorVolumeXF>
{
protected:
    size_t m_flips, m_phases, m_lines = 0;
    bool m_reorderPhase = false, m_reorderBlock = false;

public:
    typedef QI::VectorVolumeXF TInputImage;
    typedef QI::VectorVolumeXF TOutputImage;
    typedef QI::VolumeF        TMask;
    typedef MinEnergyFilter    Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;

    itkNewMacro(Self);
    itkTypeMacro(MinEnergyFilter, ImageToImageFilter);

    void setReorderPhase(const bool p) { m_reorderPhase = p; }
    void setReorderBlock(const bool b) { m_reorderBlock = b; }
    void SetPhases(const size_t p) {
        if (p < 4)
            QI_EXCEPTION("Must have a minimum of 4 phase-cycling patterns.");
        if ((p % 2) != 0)
            QI_EXCEPTION("Number of phases must be even.");
        m_phases = p;
        m_lines = m_phases / 2;
        this->Modified();
    }
    void SetInput(const TInputImage *img) ITK_OVERRIDE { this->SetNthInput(0, const_cast<TInputImage*>(img)); }
    void SetPass1(const TOutputImage *img) { this->SetNthInput(1, const_cast<TOutputImage*>(img)); }
    void SetMask(const TMask *mask) { this->SetNthInput(2, const_cast<TMask*>(mask)); }
    typename TInputImage::ConstPointer GetInput() const { return static_cast<const TInputImage *>(this->ProcessObject::GetInput(0)); }
    typename TOutputImage::ConstPointer GetPass1() const { return static_cast<const TOutputImage *>(this->ProcessObject::GetInput(1)); }
    typename TMask::ConstPointer GetMask() const { return static_cast<const TMask *>(this->ProcessObject::GetInput(2)); }

    void GenerateOutputInformation() ITK_OVERRIDE {
        Superclass::GenerateOutputInformation();
        if ((this->GetInput()->GetNumberOfComponentsPerPixel() % m_phases) != 0) {
            QI_EXCEPTION("Input size " << this->GetInput()->GetNumberOfComponentsPerPixel() << " and number of phases " << m_phases << " do not match");
        }
        m_flips = (this->GetInput()->GetNumberOfComponentsPerPixel() / m_phases);
        auto op = this->GetOutput();
        op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
        op->SetNumberOfComponentsPerPixel(m_flips);
        op->Allocate();
    }

protected:
    MinEnergyFilter() {
        this->SetNumberOfRequiredInputs(2);
        this->SetNumberOfRequiredOutputs(1);
        this->SetNthOutput(0, this->MakeOutput(0));
        this->SetPhases(4);
    }
    ~MinEnergyFilter() {}

    void DynamicThreadedGenerateData(const TInputImage::RegionType &region) ITK_OVERRIDE {
        //std::cout <<  __PRETTY_FUNCTION__ << endl;
        ConstNeighborhoodIterator<TInputImage>::RadiusType radius;
        radius.Fill(1);
        ConstNeighborhoodIterator<TInputImage> inputIter(radius, this->GetInput(), region);
        ConstNeighborhoodIterator<TInputImage> pass1Iter(radius, this->GetPass1(), region);
        auto m = this->GetMask();
        ImageRegionConstIterator<TMask> maskIter;
        if (m) {
            maskIter = ImageRegionConstIterator<TMask>(m, region);
        }
        ImageRegionIterator<TOutputImage> outputIter(this->GetOutput(), region);
        VariableLengthVector<std::complex<float>> output(m_flips);
        size_t phase_stride = m_flips;
        size_t flip_stride = 1;
        if (m_reorderPhase)
            std::swap(flip_stride, phase_stride);
        Eigen::ArrayXcd a_center(m_lines);
        Eigen::ArrayXcd b_center(m_lines);
        Eigen::ArrayXcd a_pixel(m_lines);
        Eigen::ArrayXcd b_pixel(m_lines);
        while(!inputIter.IsAtEnd()) {
            if (!m || maskIter.Get()) {
                VariableLengthVector<std::complex<float>> center_pixel = inputIter.GetCenterPixel();
                for (size_t f = 0; f < m_flips; f++) {
                    const Eigen::ArrayXcd center_array = Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>>(center_pixel.GetDataPointer() + f*flip_stride, m_phases, Eigen::InnerStride<>(phase_stride)).cast<std::complex<double>>();
                    QI::SplitBlocks(center_array, a_center, b_center, m_reorderBlock);
                    Eigen::ArrayXcd sums = Eigen::ArrayXcd::Zero(m_lines);
                    Eigen::ArrayXd nums = Eigen::ArrayXd::Zero(m_lines);
                    Eigen::ArrayXd dens = Eigen::ArrayXd::Zero(m_lines);
                    Eigen::ArrayXcd ws = Eigen::ArrayXcd::Zero(m_lines);
                    for (size_t p = 0; p < inputIter.Size(); ++p) {
                        VariableLengthVector<std::complex<float>> pass1_pixel = pass1Iter.GetPixel(p);
                        VariableLengthVector<std::complex<float>> input_pixel = inputIter.GetPixel(p);
                        const std::complex<double> Id = pass1_pixel[f];
                        if (norm(Id) > 0.) {
                            const Eigen::ArrayXcd input_array = Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>>(input_pixel.GetDataPointer() + f*flip_stride, m_phases, Eigen::InnerStride<>(phase_stride)).cast<std::complex<double>>();
                            QI::SplitBlocks(input_array, a_pixel, b_pixel, m_reorderBlock);
                            nums += real(conj(b_pixel - Id)*(b_pixel - a_pixel) + conj(b_pixel - a_pixel)*(b_pixel - Id));
                            dens += real(conj(a_pixel - b_pixel)*(a_pixel - b_pixel));
                        }
                    }
                    ws = nums / (2. * dens);
                    sums = ws*a_center + (1. - ws)*b_center;
                    output[f] = static_cast<std::complex<float>>(sums.sum() / static_cast<double>(m_lines));
                }
            } else {
                output.Fill(std::complex<float>(0.f,0.f));
            }
            outputIter.Set(output);
            ++inputIter;
            ++pass1Iter;
            ++outputIter;
            if (m) {
                ++maskIter;
            }
        }
        //std::cout << "End " << __PRETTY_FUNCTION__ );
    }

private:
    MinEnergyFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

} // End namespace itk

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();

    args::ArgumentParser parser("Removes bands from SSFP images.\n"
                                "http://github.com/spinicist/QUIT");
    
    args::Positional<std::string> input_path(parser, "INPUT", "Input filename");
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<std::string> out_arg(parser, "OUTPREFIX", "Change output prefix (default input filename)", {'o', "out"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::Flag     alt_order(parser, "ALTERNATE", "Opposing phase-incs alternate (default is 2 blocks)", {"alt-order"});
    args::Flag     ph_order(parser, "PHASE 1st", "Data order is phase, then flip-angle (default opposite)", {"ph-order"});
    args::ValueFlag<int> ph_incs(parser, "PHASE-INCS", "Number of phase increments (default 4)", {"ph-incs"}, 4);
    args::Flag     magnitude(parser, "MAGNITUDE", "Output a magnitude image only (default is complex)", {"magnitude"});
    args::ValueFlag<std::string> method(parser, "METHOD", "Choose banding-removal method. G = Geometric Solution, X = Complex Average, R = Root Mean Square, M = Maximum, N = Mean Magnitude. Default = G", {"method"},"G");
    args::ValueFlag<std::string> regularise(parser, "REGULARISE", "Chose regularisation method for GS. M = Magnitude, L = Line, N = None", {"regularise"}, "L");
    args::Flag     two_pass(parser, "SECOND PASS", "Use energy-minimisation 2nd pass scheme", {'2',"2pass"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    
    QI_LOG(verbose, "Opening input file: " << QI::CheckPos(input_path));
    auto inFile = QI::ReadVectorImage<std::complex<float>>(QI::CheckPos(input_path));
    size_t nVols = inFile->GetNumberOfComponentsPerPixel() / ph_incs.Get();
    QI_LOG(verbose, "Number of phase increments is " << ph_incs.Get() <<
                    "Number of volumes to process is " << nVols );

    std::shared_ptr<QI::BandAlgo> algo = nullptr;
    std::string suffix = "";
    if (method.Get() == "G") {
        suffix = "GS";
        QI_LOG(verbose, "Geometric solution selected" );
        auto g = std::make_shared<QI::GSAlgo>();
        g->setInputSize(inFile->GetNumberOfComponentsPerPixel());
        if (regularise.Get() == "L") {
            suffix += "L"; g->setRegularise(QI::RegEnum::Line);
        } else if (regularise.Get() == "M") {
            suffix += "M"; g->setRegularise(QI::RegEnum::Magnitude);
        } else if (regularise.Get() == "N") {
            g->setRegularise(QI::RegEnum::None);
        } else {
            QI_FAIL("Invalid regularisation " << regularise.Get() << " selected.");
        }
        algo = g;
    } else if (method.Get() == "X") {
        suffix = "CS"; algo = std::make_shared<QI::CSAlgo>();
    } else if (method.Get() == "R") {
        suffix = "RMS"; algo = std::make_shared<QI::RMSAlgo>();
    } else if (method.Get() == "N") {
        suffix = "MagMean"; algo = std::make_shared<QI::MagMeanAlgo>();
    } else if (method.Get() == "M") {
        suffix = "Max"; algo = std::make_shared<QI::MaxAlgo>();
    } else {
        QI_FAIL("Invalid method " << method.Get() << " selected.");
    }
    QI_LOG(verbose, suffix << " method selected." );
    algo->setPhases(ph_incs.Get());
    algo->setInputSize(inFile->GetNumberOfComponentsPerPixel());
    algo->setReorderPhase(ph_order);
    algo->setReorderBlock(alt_order);
    auto pass1 = QI::ApplyVectorXF::New();
    pass1->SetAlgorithm(algo);
    if (mask) pass1->SetMask(QI::ReadImage(mask.Get()));
    pass1->SetInput(0, inFile);
    pass1->SetVerbose(verbose);
    QI_LOG(verbose, "1st pass");
    if (verbose) {
        auto monitor = QI::GenericMonitor::New();
        pass1->AddObserver(itk::ProgressEvent(), monitor);
    }
    pass1->Update();
    QI::VectorVolumeXF::Pointer output = ITK_NULLPTR;
    if (two_pass) {
        suffix += "2";
        auto pass2 = itk::MinEnergyFilter::New();
        pass2->SetPhases(ph_incs.Get());
        pass2->setReorderBlock(alt_order);
        pass2->setReorderPhase(ph_order);
        pass2->SetInput(inFile);
        pass2->SetPass1(pass1->GetOutput(0));
        if (mask) pass2->SetMask(QI::ReadImage(mask.Get()));
        QI_LOG(verbose, "2nd pass");
        if (verbose) {
            auto monitor = QI::GenericMonitor::New();
            pass2->AddObserver(itk::ProgressEvent(), monitor);
        }
        pass2->Update();
        output = pass2->GetOutput();
    } else {
        output = pass1->GetOutput(0);
    }
    std::string prefix = (out_arg ? out_arg.Get() : QI::StripExt(input_path.Get()));
    std::string outname = prefix + "_" + suffix + QI::OutExt();
    QI_LOG(verbose, "Output filename: " << outname );
    if (magnitude) {
        QI::WriteVectorMagnitudeImage<QI::VectorVolumeXF>(output, outname);
    } else {
        QI::WriteVectorImage(output, outname);
    }
    QI_LOG(verbose, "Finished." );
    return EXIT_SUCCESS;
}
