/*
 *  qsignal.cpp
 *
 *  Created by Tobias Wood on 2015/06/02.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <string>
#include <iostream>
#include <Eigen/Dense>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImage.h"
#include "itkImageToImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkTimeProbe.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkProgressReporter.h"

#include "VectorToImageFilter.h"

#include "Model.h"
#include "Util.h"
#include "ThreadPool.h"
#include "Args.h"
#include "IO.h"

#include "SPGRSequence.h"
#include "SSFPSequence.h"
#include "MultiEchoSequence.h"
#include "AFISequence.h"
#include "MPRAGESequence.h"
#include "SequenceGroup.h"

/*
 * Signals Filter
 */
typedef itk::Image<float, 3> TImage;
typedef itk::VectorImage<float, 3> TVImage;
typedef itk::VectorImage<std::complex<float>, 3> TCVImage;

class SignalsFilter : public itk::ImageToImageFilter<TImage, TCVImage> {
protected:
    std::shared_ptr<QI::SequenceBase> m_sequence;
    std::shared_ptr<QI::Model> m_model;
    double m_sigma = 0.0;
    itk::TimeProbe m_clock;
    itk::RealTimeClock::TimeStampType m_meanTime = 0.0, m_totalTime = 0.0;
    itk::SizeValueType m_evaluations = 0;

public:
    /** Standard class typedefs. */
    typedef SignalsFilter                        Self;
    typedef ImageToImageFilter<TImage, TCVImage> Superclass;
    typedef itk::SmartPointer<Self>              Pointer;
    typedef typename TImage::RegionType          TRegion;

    itkNewMacro(Self); /** Method for creation through the object factory. */
    itkTypeMacro(Self, Superclass); /** Run-time type information (and related methods). */

    void SetInput(const size_t i, const TImage *img) {
        if (i < m_model->nParameters()) {
            this->SetNthInput(i, const_cast<TImage*>(img));
        } else {
            QI_EXCEPTION("Const input " << i << " out of range");
        }
    }
    void SetMask(const TImage *mask) {
        this->SetNthInput(m_model->nParameters(), const_cast<TImage*>(mask));
    }
    
    typename TImage::ConstPointer GetInput(const size_t i) const {
        if (i < m_model->nParameters()) {
            return static_cast<const TImage *> (this->ProcessObject::GetInput(i));
        } else {
            QI_EXCEPTION("Get Data Input " << i << " out of range.");
        }
    }

    typename TImage::ConstPointer GetMask() const {
        return static_cast<const TImage *>(this->ProcessObject::GetInput(m_model->nParameters()));
    }

    TCVImage *GetOutput() {
        return dynamic_cast<TCVImage *>(this->ProcessObject::GetOutput(0));
    }

    void SetSequence(std::shared_ptr<QI::SequenceBase> s) {
        m_sequence = s;
        this->SetNumberOfRequiredOutputs(1);
        this->SetNthOutput(0, this->MakeOutput(0));
    }
    void SetModel(std::shared_ptr<QI::Model> m) {
        m_model = m;
        this->SetNumberOfRequiredInputs(1);
    }
    void SetSigma(const double s) { m_sigma = s; }

    itk::RealTimeClock::TimeStampType GetTotalTime() const { return m_totalTime; }
    itk::RealTimeClock::TimeStampType GetMeanTime() const { return m_meanTime; }
    itk::SizeValueType GetEvaluations() const { return m_evaluations; }

    void GenerateOutputInformation() ITK_OVERRIDE {
        //std::cout <<  __PRETTY_FUNCTION__ << std::endl;
        Superclass::GenerateOutputInformation();
        const auto op = this->GetOutput();
        const auto ip = this->GetInput(0);
        op->SetRegions(ip->GetLargestPossibleRegion());
        op->SetSpacing(ip->GetSpacing());
        op->SetOrigin(ip->GetOrigin());
        op->SetDirection(ip->GetDirection());
        op->SetNumberOfComponentsPerPixel(m_sequence->size());
        op->Allocate(true);
    }

protected:
    SignalsFilter() {}
    ~SignalsFilter(){}

    void ThreadedGenerateData(const OutputImageRegionType &region, itk::ThreadIdType threadId) ITK_OVERRIDE {
        std::vector<itk::ImageRegionConstIterator<TImage>> inIters(m_model->nParameters());
        for (size_t i = 0; i < m_model->nParameters(); i++) {
            if (this->GetInput(i))
                inIters[i] = itk::ImageRegionConstIterator<TImage>(this->GetInput(i), region);
        }
        typename TImage::ConstPointer mask = this->GetMask();
        itk::ImageRegionConstIterator<TImage> maskIter;
        if (mask) {
            maskIter = itk::ImageRegionConstIterator<TImage>(mask, region);
        }
        itk::ImageRegionIterator<TCVImage> outputIter(this->GetOutput(), region);
        
        if (threadId == 0)
            m_clock.Reset();
        while(!outputIter.IsAtEnd()) {
            if (!mask || maskIter.Get()) {
                if (threadId == 0) {
                    m_clock.Start();
                }
                Eigen::VectorXd parameters = m_model->Default();
                for (size_t i = 0; i < inIters.size(); i++) {
                    if (this->GetInput(i))
                        parameters[i] = inIters[i].Get();
                }
                Eigen::VectorXcd allData = m_sequence->signal(m_model, parameters);
                if (m_sigma != 0.0) {
                    Eigen::VectorXcd noise(m_sequence->size());
                    // Simple Box Muller transform
                    Eigen::ArrayXd U = (Eigen::ArrayXd::Random(m_sequence->size()) * 0.5) + 0.5;
                    Eigen::ArrayXd V = (Eigen::ArrayXd::Random(m_sequence->size()) * 0.5) + 0.5;
                    noise.real() = (m_sigma / M_SQRT2) * (-2. * U.log()).sqrt() * cos(2. * M_PI * V);
                    noise.imag() = (m_sigma / M_SQRT2) * (-2. * V.log()).sqrt() * sin(2. * M_PI * U);
                    allData += noise;
                }
                Eigen::VectorXcf floatData = allData.cast<std::complex<float>>();
                itk::VariableLengthVector<std::complex<float>> dataVector(floatData.data(), m_sequence->size());
                outputIter.Set(dataVector);
                if (threadId == 0) {
                    m_clock.Stop();
                }
            }
            if (mask) {
                ++maskIter;
            }
            for (size_t i = 0; i < m_model->nParameters(); i++) {
                if (this->GetInput(i))
                    ++inIters[i];
            }
            ++outputIter;
        }
        if (threadId == 0) {
            m_evaluations = m_clock.GetNumberOfStops();
            m_meanTime = m_clock.GetMean();
            m_totalTime = m_clock.GetTotal();
        }
    }

    itk::DataObject::Pointer MakeOutput(unsigned int idx) {
        itk::DataObject::Pointer output;
        if (idx < 1) {
            auto img = TCVImage::New();
            output = img;
        } else {
            std::cerr << "No output " << idx << std::endl;
            output = NULL;
        }
        return output.GetPointer();
    }

private:
    SignalsFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Generates simulated images from signal equations.\n"
                                "Input is read from stdin as for analysis programs.\n"
                                "http://github.com/spinicist/QUIT");
    args::PositionalList<std::string> filenames(parser, "OUTPUT FILES", "Output file names");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> ref_arg(parser, "REFERENCE", "Resample inputs to this reference space", {'r', "ref"});
    args::ValueFlag<float> noise(parser, "NOISE", "Add complex noise with std=value", {'N',"noise"}, 0.);
    args::ValueFlag<int> seed(parser, "SEED", "Seed noise RNG with specific value", {'s', "seed"}, -1);
    args::ValueFlag<int> model_arg(parser, "MODEL", "Choose number of components in model (1/2/3)", {'M',"model"}, 1);
    args::Flag     complex(parser, "COMPLEX", "Save complex images", {'x',"complex"});
    QI::ParseArgs(parser, argc, argv);

    std::shared_ptr<QI::Model> model = nullptr;
    switch (model_arg.Get()) {
        case 1: model = std::make_shared<QI::SCD>(); break;
        case 2: model = std::make_shared<QI::MCD2>(); break;
        case 3: model = std::make_shared<QI::MCD3>(); break;
        default:
            QI_FAIL("Unknown number of components: " << model_arg.Get());
    }
    if (verbose) std::cout << "Using " << model->Name() << " model." << std::endl;

    if (seed) {
        std::srand(seed.Get());
    } else {
        std::srand((unsigned int) time(0));
    }

    /***************************************************************************
     * Read in parameter files
     **************************************************************************/
    SignalsFilter::Pointer calcSignal = SignalsFilter::New();
    calcSignal->SetModel(model);
    calcSignal->SetSigma(noise.Get());
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(threads.Get());
    if (mask) {
        if (verbose) std::cout << "Reading mask: " << mask.Get() << std::endl;
        calcSignal->SetMask(QI::ReadImage(mask.Get()));
    }

    QI::VolumeF::Pointer reference = ITK_NULLPTR;
    if (ref_arg) {
        if (verbose) std::cout << "Reading reference image: " << ref_arg.Get() << std::endl;
        reference = QI::ReadImage(ref_arg.Get());
    }
    if (verbose) {
        auto monitor = QI::GenericMonitor::New();
        calcSignal->AddObserver(itk::ProgressEvent(), monitor);
    }
    if (verbose) std::cout << "Loading parameters." << std::endl;
    cereal::JSONInputArchive input(std::cin);
    for (size_t i = 0; i < model->nParameters(); i++) {
        std::string par_filename;
        std::string par_name = model->ParameterNames()[i];
        input(cereal::make_nvp(par_name, par_filename));
        if (par_filename != "") {
            if (verbose) std::cout << "Opening " << par_filename << std::endl;
            QI::VolumeF::Pointer param = QI::ReadImage(par_filename);
            if (reference) {
                if (verbose) std::cout << "Resampling to reference" << std::endl;
                typedef itk::ResampleImageFilter<QI::VolumeF, QI::VolumeF, double> TResampler;
                typedef itk::LinearInterpolateImageFunction<QI::VolumeF, double> TInterp;
                typename TInterp::Pointer interp = TInterp::New();
                interp->SetInputImage(param);
                typename TResampler::Pointer resamp = TResampler::New();
                resamp->SetInput(param);
                resamp->SetInterpolator(interp);
                resamp->SetDefaultPixelValue(0.);
                resamp->SetOutputParametersFromImage(reference);
                resamp->Update();
                QI::VolumeF::Pointer rparam = resamp->GetOutput();
                rparam->DisconnectPipeline();   
                calcSignal->SetInput(i, rparam);
            } else { 
                calcSignal->SetInput(i, param);
            }
        } else {
            if (verbose) std::cout << "Using default " << par_name << " value: " << model->Default()[i] << std::endl;
        }
    }

    /***************************************************************************
     * Set up sequences
     **************************************************************************/
    auto sequences = QI::ReadSequence<QI::SequenceGroup>(input, false);
    if (filenames.Get().size() != sequences.count()) {
        QI_FAIL("Input filenames size " << filenames.Get().size()
                << " does not match sequences size " << sequences.count());
    } else {
        if (verbose) std::cout << "Found " << sequences.count() << " sequences to generate" << std::endl;
    }
    for (size_t i = 0; i < sequences.count(); i++) {
        if (verbose) {
            std::cout << "Sequence: " << std::endl;
            cereal::JSONOutputArchive archive(std::cout);
            archive(sequences[i]);
        }
        calcSignal->SetSequence(sequences[i]);
        calcSignal->Update();
        QI::VectorVolumeXF::Pointer output = calcSignal->GetOutput();
        output->DisconnectPipeline();
        if (verbose) std::cout << "Mean evaluation time: " << calcSignal->GetMeanTime() << " s ( " << calcSignal->GetEvaluations() << " voxels)" << std::endl;
        if (complex) {
            if (verbose) std::cout << "Saving complex image: " << filenames.Get()[i] << std::endl;
            QI::WriteVectorImage(output, outarg.Get() + filenames.Get()[i]);
        } else {
            if (verbose) std::cout << "Saving magnitude image: " << filenames.Get()[i] << std::endl;
            QI::WriteVectorMagnitudeImage(output, outarg.Get() + filenames.Get()[i]);
        }
    }
    if (verbose) std::cout << "Finished all sequences." << std::endl;
    return EXIT_SUCCESS;
}

