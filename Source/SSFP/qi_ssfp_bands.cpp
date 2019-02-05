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

#include "itkComplexToModulusImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkUnaryFunctorImageFilter.h"

#include "Args.h"
#include "Banding.h"
#include "ImageIO.h"
#include "Log.h"
#include "Monitor.h"
#include "Util.h"

namespace itk {

class MinEnergyFilter : public ImageToImageFilter<QI::VectorVolumeXF, QI::VectorVolumeXF> {
  protected:
    size_t m_flips, m_phases, m_lines = 0;
    bool   m_reorderPhase = false, m_reorderBlock = false;

  public:
    typedef QI::VectorVolumeXF                            TInputImage;
    typedef QI::VectorVolumeXF                            TOutputImage;
    typedef QI::VolumeF                                   TMask;
    typedef MinEnergyFilter                               Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self>                            Pointer;

    itkNewMacro(Self);
    itkTypeMacro(MinEnergyFilter, ImageToImageFilter);

    void setReorderPhase(const bool p) { m_reorderPhase = p; }
    void setReorderBlock(const bool b) { m_reorderBlock = b; }
    void SetPhases(const size_t p) {
        if (p < 4)
            QI::Fail("Must have a minimum of 4 phase-cycling patterns.");
        if ((p % 2) != 0)
            QI::Fail("Number of phases must be even.");
        m_phases = p;
        m_lines  = m_phases / 2;
        this->Modified();
    }
    void SetInput(const TInputImage *img) ITK_OVERRIDE {
        this->SetNthInput(0, const_cast<TInputImage *>(img));
    }
    void SetPass1(const TOutputImage *img) {
        this->SetNthInput(1, const_cast<TOutputImage *>(img));
    }
    void SetMask(const TMask *mask) { this->SetNthInput(2, const_cast<TMask *>(mask)); }
    typename TInputImage::ConstPointer GetInput() const {
        return static_cast<const TInputImage *>(this->ProcessObject::GetInput(0));
    }
    typename TOutputImage::ConstPointer GetPass1() const {
        return static_cast<const TOutputImage *>(this->ProcessObject::GetInput(1));
    }
    typename TMask::ConstPointer GetMask() const {
        return static_cast<const TMask *>(this->ProcessObject::GetInput(2));
    }

    void GenerateOutputInformation() ITK_OVERRIDE {
        Superclass::GenerateOutputInformation();
        if ((this->GetInput()->GetNumberOfComponentsPerPixel() % m_phases) != 0) {
            QI::Fail("Input size {} and number of phase {} do not match",
                     this->GetInput()->GetNumberOfComponentsPerPixel(), m_phases);
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
        // std::cout <<  __PRETTY_FUNCTION__ << endl;
        ConstNeighborhoodIterator<TInputImage>::RadiusType radius;
        radius.Fill(1);
        ConstNeighborhoodIterator<TInputImage> inputIter(radius, this->GetInput(), region);
        ConstNeighborhoodIterator<TInputImage> pass1Iter(radius, this->GetPass1(), region);
        auto                                   m = this->GetMask();
        ImageRegionConstIterator<TMask>        maskIter;
        if (m) {
            maskIter = ImageRegionConstIterator<TMask>(m, region);
        }
        ImageRegionIterator<TOutputImage>         outputIter(this->GetOutput(), region);
        VariableLengthVector<std::complex<float>> output(m_flips);
        size_t                                    phase_stride = m_flips;
        size_t                                    flip_stride  = 1;
        if (m_reorderPhase)
            std::swap(flip_stride, phase_stride);
        Eigen::ArrayXcd a_center(m_lines);
        Eigen::ArrayXcd b_center(m_lines);
        Eigen::ArrayXcd a_pixel(m_lines);
        Eigen::ArrayXcd b_pixel(m_lines);
        while (!inputIter.IsAtEnd()) {
            if (!m || maskIter.Get()) {
                VariableLengthVector<std::complex<float>> center_pixel = inputIter.GetCenterPixel();
                for (size_t f = 0; f < m_flips; f++) {
                    const Eigen::ArrayXcd center_array =
                        Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>>(
                            center_pixel.GetDataPointer() + f * flip_stride, m_phases,
                            Eigen::InnerStride<>(phase_stride))
                            .cast<std::complex<double>>();
                    QI::SplitBlocks(center_array, a_center, b_center, m_reorderBlock);
                    Eigen::ArrayXcd sums = Eigen::ArrayXcd::Zero(m_lines);
                    Eigen::ArrayXd  nums = Eigen::ArrayXd::Zero(m_lines);
                    Eigen::ArrayXd  dens = Eigen::ArrayXd::Zero(m_lines);
                    Eigen::ArrayXcd ws   = Eigen::ArrayXcd::Zero(m_lines);
                    for (size_t p = 0; p < inputIter.Size(); ++p) {
                        VariableLengthVector<std::complex<float>> pass1_pixel =
                            pass1Iter.GetPixel(p);
                        VariableLengthVector<std::complex<float>> input_pixel =
                            inputIter.GetPixel(p);
                        const std::complex<double> Id = pass1_pixel[f];
                        if (norm(Id) > 0.) {
                            const Eigen::ArrayXcd input_array =
                                Eigen::Map<const Eigen::ArrayXcf, 0, Eigen::InnerStride<>>(
                                    input_pixel.GetDataPointer() + f * flip_stride, m_phases,
                                    Eigen::InnerStride<>(phase_stride))
                                    .cast<std::complex<double>>();
                            QI::SplitBlocks(input_array, a_pixel, b_pixel, m_reorderBlock);
                            nums += real(conj(b_pixel - Id) * (b_pixel - a_pixel) +
                                         conj(b_pixel - a_pixel) * (b_pixel - Id));
                            dens += real(conj(a_pixel - b_pixel) * (a_pixel - b_pixel));
                        }
                    }
                    ws   = nums / (2. * dens);
                    sums = ws * a_center + (1. - ws) * b_center;
                    output[f] =
                        static_cast<std::complex<float>>(sums.sum() / static_cast<double>(m_lines));
                }
            } else {
                output.Fill(std::complex<float>(0.f, 0.f));
            }
            outputIter.Set(output);
            ++inputIter;
            ++pass1Iter;
            ++outputIter;
            if (m) {
                ++maskIter;
            }
        }
        // std::cout << "End " << __PRETTY_FUNCTION__ );
    }

  private:
    MinEnergyFilter(const Self &); // purposely not implemented
    void operator=(const Self &);  // purposely not implemented
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
    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<std::string> out_arg(
        parser, "OUTPREFIX", "Change output prefix (default input filename)", {'o', "out"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"}, 4);
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask",
                                      {'m', "mask"});
    args::Flag alt_order(parser, "ALTERNATE", "Opposing phase-incs alternate (default is 2 blocks)",
                         {"alt-order"});
    args::Flag ph_order(parser, "PHASE 1st",
                        "Data order is phase, then flip-angle (default opposite)", {"ph-order"});
    args::ValueFlag<int> ph_incs(parser, "PHASE-INCS", "Number of phase increments (default 4)",
                                 {"ph-incs"}, 4);
    args::Flag magnitude(parser, "MAGNITUDE", "Output a magnitude image only (default is complex)",
                         {"magnitude"});
    args::ValueFlag<std::string> method(
        parser, "METHOD",
        "Choose banding-removal method. G = Geometric Solution, X = Complex Average, R = Root Mean "
        "Square, M = Maximum, N = Mean Magnitude. Default = G",
        {"method"}, "G");
    args::ValueFlag<std::string> regularise(
        parser, "REGULARISE",
        "Chose regularisation method for GS. M = Magnitude, L = Line, N = None", {"regularise"},
        "L");
    args::Flag two_pass(parser, "SECOND PASS", "Use energy-minimisation 2nd pass scheme",
                        {'2', "2pass"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI::Log(verbose, "Opening input file: {}", QI::CheckPos(input_path));
    auto         inFile = QI::ReadVectorImage<std::complex<float>>(QI::CheckPos(input_path));
    const size_t nVols  = inFile->GetNumberOfComponentsPerPixel();
    QI::Log(verbose, "Phase increments = {} Number of volumes = {}", ph_incs.Get(), nVols);
    QI::VectorVolumeXF::Pointer output = ITK_NULLPTR;
    std::string                 suffix = "";
    if (method.Get() == "G") {
        suffix = "GS";
        QI::Log(verbose, "Geometric solution selected");
        QI::GSFunctor gs(nVols, ph_incs.Get(), ph_order, alt_order);
        if (regularise.Get() == "L") {
            suffix += "L";
            gs.setRegularise(QI::RegEnum::Line);
        } else if (regularise.Get() == "M") {
            suffix += "M";
            gs.setRegularise(QI::RegEnum::Magnitude);
        } else if (regularise.Get() == "N") {
            gs.setRegularise(QI::RegEnum::None);
        } else {
            QI::Fail("Invalid regularisation {}", regularise.Get());
        }
        QI::Log(verbose, "Regularisation = {}", suffix);
        auto pass1 = itk::UnaryFunctorImageFilter<QI::VectorVolumeXF, QI::VectorVolumeXF,
                                                  QI::GSFunctor>::New();
        pass1->SetFunctor(gs);
        pass1->SetInput(inFile);
        pass1->Update();
        output = pass1->GetOutput();
    } else if (method.Get() == "X") {
        suffix = "CS";
        QI::CSFunctor functor(nVols, ph_incs.Get(), ph_order, alt_order);
        auto          pass1 = itk::UnaryFunctorImageFilter<QI::VectorVolumeXF, QI::VectorVolumeXF,
                                                  QI::CSFunctor>::New();
        pass1->SetFunctor(functor);
        pass1->SetInput(inFile);
        pass1->Update();
        output = pass1->GetOutput();
    } else if (method.Get() == "R") {
        suffix = "RMS";
        QI::RMSFunctor functor(nVols, ph_incs.Get(), ph_order, alt_order);
        auto           pass1 = itk::UnaryFunctorImageFilter<QI::VectorVolumeXF, QI::VectorVolumeXF,
                                                  QI::RMSFunctor>::New();
        pass1->SetFunctor(functor);
        pass1->SetInput(inFile);
        pass1->Update();
        output = pass1->GetOutput();
    } else if (method.Get() == "N") {
        suffix = "MagMean";
        QI::MagMeanFunctor functor(nVols, ph_incs.Get(), ph_order, alt_order);
        auto pass1 = itk::UnaryFunctorImageFilter<QI::VectorVolumeXF, QI::VectorVolumeXF,
                                                  QI::MagMeanFunctor>::New();
        pass1->SetFunctor(functor);
        pass1->SetInput(inFile);
        pass1->Update();
        output = pass1->GetOutput();
    } else if (method.Get() == "M") {
        suffix = "Max";
        QI::MaxFunctor functor(nVols, ph_incs.Get(), ph_order, alt_order);
        auto           pass1 = itk::UnaryFunctorImageFilter<QI::VectorVolumeXF, QI::VectorVolumeXF,
                                                  QI::MaxFunctor>::New();
        pass1->SetFunctor(functor);
        pass1->SetInput(inFile);
        pass1->Update();
        output = pass1->GetOutput();
    } else {
        QI::Fail("Invalid method: {}", method.Get());
    }

    if (two_pass) {
        suffix += "2";
        auto pass2 = itk::MinEnergyFilter::New();
        pass2->SetPhases(ph_incs.Get());
        pass2->setReorderBlock(alt_order);
        pass2->setReorderPhase(ph_order);
        pass2->SetInput(inFile);
        pass2->SetPass1(output);
        if (mask)
            pass2->SetMask(QI::ReadImage(mask.Get()));
        QI::Log(verbose, "2nd pass");
        if (verbose) {
            auto monitor = QI::GenericMonitor::New();
            pass2->AddObserver(itk::ProgressEvent(), monitor);
        }
        pass2->Update();
        output = pass2->GetOutput();
    }
    std::string prefix  = (out_arg ? out_arg.Get() : QI::StripExt(input_path.Get()));
    std::string outname = prefix + "_" + suffix + QI::OutExt();
    QI::Log(verbose, "Writing output file: {}", outname);
    if (magnitude) {
        QI::WriteVectorMagnitudeImage<QI::VectorVolumeXF>(output, outname);
    } else {
        QI::WriteVectorImage(output, outname);
    }
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
