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
#include "Eigen/Dense"

#include "itkUnaryFunctorImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkConstNeighborhoodIterator.h"

#include "Filters/ReorderVectorFilter.h"
#include "QI/Util.h"
#include "QI/Filters/GeometricSolutionFilter.h"
#include "QI/Option.h"

using namespace std;
using namespace Eigen;

class BandAlgo : public QI::ApplyVectorXF::Algorithm {
protected:
    size_t m_flips, m_lines, m_crossings, m_phases = 4;
    bool m_reorderPhase = false, m_reorderBlock = false;
    TOutput m_zero;
public:
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 0; }
    size_t numOutputs() const override { return 1; }
    size_t dataSize() const override { return m_flips * m_phases; }
    size_t outputSize(const int i) const override { return m_flips; }
    void setPhases(const size_t p) {
        if (p < 4)
            QI_EXCEPTION("Must have a minimum of 4 phase-cycling patterns.");
        if ((p % 2) != 0)
            QI_EXCEPTION("Number of phases must be even.");
        m_phases = p;
        m_lines = m_phases / 2;
        m_crossings = QI::Choose(m_lines, 2);
    }
    void setInputSize(const size_t s) {
        m_flips = s / m_phases;
        m_zero.SetSize(m_flips);
        m_zero.Fill(complex<float>(0.));
    }
    void setReorderPhase(const bool p) { m_reorderPhase = p; }
    void setReorderBlock(const bool b) { m_reorderBlock = b; }

    virtual std::vector<float> defaultConsts() override {
        std::vector<float> def;
        return def;
    }
    virtual const TOutput &zero(const size_t i) const override {
        return m_zero;
    }
    virtual complex<float> applyFlip(const Map<const ArrayXcf, 0, InnerStride<>> &vf) const = 0;
    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                       std::vector<TOutput> &outputs, TConst &residual,
                       TInput &resids, TIters &its) const override
    {
        size_t phase_stride = 1;
        if (m_reorderPhase)
            phase_stride = m_flips;
        ArrayXcf out;
        for (int f = 0; f < m_flips; f++) {
            const Map<const ArrayXcf, 0, InnerStride<>> vf(inputs[0].GetDataPointer() + f, inputs[0].GetSize(), InnerStride<>(phase_stride));
            outputs[0][f] = this->applyFlip(vf);
        }
    }
};

class CSAlgo : public BandAlgo {
public:
    complex<float> applyFlip(const Map<const ArrayXcf, 0, InnerStride<>> &vf) const override {
        return vf.mean();
    }
};

class MagMeanAlgo : public BandAlgo {
public:
    complex<float> applyFlip(const Map<const ArrayXcf, 0, InnerStride<>> &vf) const override {
        return vf.abs().mean();
    }
};

class RMSAlgo : public BandAlgo {
public:
    complex<float> applyFlip(const Map<const ArrayXcf, 0, InnerStride<>> &vf) const override {
        float sum = vf.abs().square().sum();
        return complex<float>(sqrt(sum / vf.rows()), 0.);
    }
};

class MaxAlgo : public BandAlgo {
public:
    complex<float> applyFlip(const Map<const ArrayXcf, 0, InnerStride<>> &vf) const override {
        complex<float> max = std::numeric_limits<complex<float>>::lowest();
        for (size_t i = 0; i < vf.rows(); i++) {
            if (std::abs(vf[i]) > std::abs(max)) max = vf[i];
        }
        return max;
    }
};

class GSAlgo : public BandAlgo {
public:
    enum class RegEnum { None = 0, Line, Magnitude };
protected:
    RegEnum m_Regularise = RegEnum::Line;
public:
    const RegEnum &regularise()       { return m_Regularise; }
    void setRegularise(const RegEnum &r) { m_Regularise = r;}
    complex<float> applyFlip(const Map<const ArrayXcf, 0, InnerStride<>> &vf) const override {
        Eigen::ArrayXcd a(m_flips / 2);
        Eigen::ArrayXcd b(m_flips / 2);
        if (m_reorderBlock) {
            for (int i = 0; i < m_phases; i++) {
                a[i] = static_cast<complex<double>>(vf[i*2]);
                b[i] = static_cast<complex<double>>(vf[i*2+1]);
            }
        } else {
            a = vf.head(m_lines).cast<std::complex<double>>();
            b = vf.tail(m_lines).cast<std::complex<double>>();
        }

        std::complex<double> sum(0., 0.);
        for (size_t i = 0; i < m_lines; i++) {
            for (size_t j = i + 1; j < m_lines; j++) {
                const std::complex<double> di = b[i] -  a[i], dj = b[j] - a[j];
                const std::complex<double> ni(di.imag(), -di.real()), nj(dj.imag(), -dj.real());

                const double mu = QI::cdot(a[j] - a[i], nj) / QI::cdot(di, nj);
                const double nu = QI::cdot(a[i] - a[j], ni) / QI::cdot(dj, ni);
                const double xi = 1.0 - pow(QI::cdot(di, dj) / (abs(di)*abs(dj)), 2.0);

                const std::complex<double> cs = (a[i] + a[j] + b[i] + b[j]) / 4.0;
                const std::complex<double> gs = a[i] + mu * di;

                switch (m_Regularise) {
                case RegEnum::None: sum += gs; break;
                case RegEnum::Magnitude:
                    if (norm(gs) < std::max({std::norm(a[i]), std::norm(a[j]), std::norm(b[i]), std::norm(b[j])})) {
                        sum += gs;
                    } else {
                        sum += cs;
                    }
                    break;
                case RegEnum::Line:
                    if ((mu > -xi) && (mu < 1 + xi) && (nu > -xi) && (nu < 1 + xi)) {
                        sum += gs;
                    } else {
                        sum += cs;
                    }
                    break;
                }
            }
        }
        return static_cast<std::complex<float>>(sum / static_cast<double>(m_crossings));
    }
};

namespace itk {

class MinEnergyFilter : public ImageToImageFilter<VectorImage<complex<float>, 3>, Image<complex<float>, 3>>
{
protected:
    size_t m_flips, m_phases, m_lines = 0;

public:
    typedef VectorImage<complex<float>, 3>     TInputImage;
    typedef Image<complex<float>, 3>           TOutputImage;
    typedef Image<float, 3>                    TMask;
    typedef MinEnergyFilter                   Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self>                 Pointer;

    itkNewMacro(Self);
    itkTypeMacro(MinEnergyFilter, ImageToImageFilter);

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

    virtual void GenerateOutputInformation() ITK_OVERRIDE {
        Superclass::GenerateOutputInformation();
        if ((this->GetInput()->GetNumberOfComponentsPerPixel() % m_phases) != 0) {
            QI_EXCEPTION("Input size and number of phases do not match");
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

    virtual void ThreadedGenerateData(const TInputImage::RegionType &region, ThreadIdType threadId) ITK_OVERRIDE {
        //std::cout <<  __PRETTY_FUNCTION__ << endl;
        ConstNeighborhoodIterator<TInputImage>::RadiusType radius;
        radius.Fill(1);
        ConstNeighborhoodIterator<TInputImage> inputIter(radius, this->GetInput(), region);
        ConstNeighborhoodIterator<TOutputImage> pass1Iter(radius, this->GetPass1(), region);

        auto m = this->GetMask();
        ImageRegionConstIterator<TMask> maskIter;
        if (m) {
            maskIter = ImageRegionConstIterator<TMask>(m, region);
        }
        ImageRegionIterator<TOutputImage> outputIter(this->GetOutput(), region);
        while(!inputIter.IsAtEnd()) {
            if (!m || maskIter.Get()) {

                complex<double> output = 0.;
                VariableLengthVector<complex<float>> inputVector = inputIter.GetCenterPixel();
                Map<const ArrayXcf> allInput(inputVector.GetDataPointer(), m_phases);
                ArrayXcd a_center = allInput.head(m_lines).cast<complex<double>>();
                ArrayXcd b_center = allInput.tail(m_lines).cast<complex<double>>();
                for (int i = 0; i < m_lines; i++) {
                    double num = 0., den = 0.;
                    for (int p = 0; p < inputIter.Size(); ++p) {
                        const complex<double> Id = static_cast<complex<double>>(pass1Iter.GetPixel(p));
                        VariableLengthVector<complex<float>> inputVector = inputIter.GetPixel(p);
                        Map<const VectorXcf> allInput(inputVector.GetDataPointer(), m_phases);
                        const ArrayXcd a_pixel = allInput.head(m_lines).cast<complex<double>>();
                        const ArrayXcd b_pixel = allInput.tail(m_lines).cast<complex<double>>();

                        num += real(conj(b_pixel[i] - Id)*(b_pixel[i] - a_pixel[i]) + conj(b_pixel[i] - a_pixel[i])*(b_pixel[i] - Id));
                        den += real(conj(a_pixel[i] - b_pixel[i])*(a_pixel[i] - b_pixel[i]));
                    }
                    double w = num / (2. * den);
                    if (isfinite(w)) {
                        output += w*a_center[i] + (1. - w)*b_center[i];
                    }
                }
                outputIter.Set(static_cast<complex<float>>(output / static_cast<double>(m_lines)));
            }
            ++inputIter;
            ++pass1Iter;
            if (m)
                ++maskIter;
            ++outputIter;
        }
        //std::cout << "End " << __PRETTY_FUNCTION__ << std::endl;
    }

private:
    MinEnergyFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

} // End namespace itk

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::OptionList  opts("Usage is: qissfpbands [options] input\n\nInput must be a single complex image with >2 pairs phase incs\n");
    QI::Option<int> num_threads(4,'T',"threads","Use N threads (default=4, 0=hardware limit)", opts);
    QI::Option<int> ph_incs(4,'\0',"ph_incs","Number of phase increments (default is 4).", opts);
    QI::Switch      ph_order('\0',"ph_order","Data order is phase, then flip-angle (default opposite).", opts);
    QI::Switch      two_pass('2',"2pass","Use the energy-minimisation scheme from Xiang and Hoff.", opts);
    QI::EnumOption  regularise("MLN",'L','R',"regularise","Apply regularisation (magnitude/line/none)", opts);
    QI::EnumOption  method("GXRMN",'G','M',"method","Choose method GXRMN = GS/CS/RMS/Max/Mag Mean", opts);
    QI::Switch      magnitude('M',"magnitude","Output a magnitude image (default is complex)", opts);
    QI::Switch      alt_order('\0',"alt_order","Opposing phase-incs alternate (default is 2 blocks)", opts);
    QI::ImageOption<QI::VolumeF> mask('m', "mask", "Mask input with specified file", opts);
    QI::Option<string> prefix("",'o',"out","Specify output filename (default input+_nobands)", opts);
    QI::Switch      verbose('v',"verbose","Print more information", opts);
    QI::Help        help(opts);
    std::vector<std::string> nonopts = opts.parse(argc, argv);
    if (nonopts.size() != 1) {
        cerr << opts << endl;
        cerr << "Incorrect number of arguments." << endl;
        return EXIT_FAILURE;
    }
    if (*verbose) cout << "Opening input file: " << nonopts[0] << endl;
    auto inFile = QI::ReadVectorImage<complex<float>>(nonopts[0]);
    size_t nVols = inFile->GetNumberOfComponentsPerPixel() / *ph_incs;
    if (*verbose) {
        cout << "Number of phase increments is " << *ph_incs << endl;
        cout << "Number of volumes to process is " << nVols << endl;
    }

    shared_ptr<BandAlgo> algo = nullptr;
    string suffix = "";
    switch (*method) {
        case 'G': {
            suffix = "GS";
            if (*verbose) cout << "Geometric solution selected" << endl;
            auto g = make_shared<GSAlgo>();
            g->setInputSize(inFile->GetNumberOfComponentsPerPixel());
            g->setReorderBlock(*alt_order);
            switch(*regularise) {
                case 'L': g->setRegularise(GSAlgo::RegEnum::Line); break;
                case 'M': g->setRegularise(GSAlgo::RegEnum::Magnitude); break;
                case 'N': g->setRegularise(GSAlgo::RegEnum::None); break;
            }
            algo = g;
        }   break;
        case 'X': suffix = "CS"; algo = make_shared<CSAlgo>(); break;
        case 'R': suffix = "RMS"; algo = make_shared<RMSAlgo>(); break;
        case 'N': suffix = "MagMean"; algo = make_shared<MagMeanAlgo>(); break;
        case 'M': suffix = "Max"; algo = make_shared<MaxAlgo>(); break;
    }
    if (*verbose) cout << suffix << " method selected." << endl;
    algo->setInputSize(inFile->GetNumberOfComponentsPerPixel());
    algo->setPhases(*ph_incs);
    algo->setReorderPhase(*ph_order);
    auto apply = QI::ApplyVectorXF::New();
    apply->SetAlgorithm(algo);
    apply->SetMask(*mask);
    apply->SetInput(0, inFile);
    apply->SetPoolsize(*num_threads);
    if (*verbose) {
        cout << "Processing" << endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    /*    if (do_2pass) {
            auto p2 = itk::MinEnergyFilter::New();
            p2->SetInput(blockVector->GetOutput());
            p2->SetPass1(gs->GetOutput());
            if (mask)
                p2->SetMask(mask);
            process = p2;
        } else {
            process = gs;
        }*/
    if (*prefix == "")
        *prefix = QI::StripExt(nonopts[0]) + "_" + suffix;
    string outname = *prefix;
    outname.append(QI::OutExt());
    if (*verbose) cout << "Output filename: " << outname << endl;
    cout << apply->GetOutput(0)->GetDirection() << endl;
    if (*magnitude) {
        QI::WriteVectorMagnitudeImage<QI::VectorVolumeXF>(apply->GetOutput(0), outname);
    } else {
        QI::WriteVectorImage(apply->GetOutput(0), outname);
    }
    if (*verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
