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
#include "itkClampImageFilter.h"

#include "Filters/VectorToImageFilter.h"

#include "QI/Models/Model.h"
#include "QI/Sequences/Sequences.h"
#include "QI/Util.h"
#include "QI/ThreadPool.h"
#include "QI/Option.h"

using namespace std;
using namespace Eigen;
using namespace QI;

//******************************************************************************
// Filter
//******************************************************************************
typedef itk::Image<float, 3> TImage;
typedef itk::VectorImage<float, 3> TVImage;
typedef itk::VectorImage<complex<float>, 3> TCVImage;

class SignalsFilter : public itk::ImageToImageFilter<TImage, TCVImage> {
protected:
    shared_ptr<SequenceBase> m_sequence;
    shared_ptr<Model> m_model;
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

    void SetSequence(shared_ptr<SequenceBase> s) {
        m_sequence = s;
        this->SetNumberOfRequiredOutputs(1);
        this->SetNthOutput(0, this->MakeOutput(0));
    }
    void SetModel(shared_ptr<Model> m) {
        m_model = m;
        this->SetNumberOfRequiredInputs(1);
    }
    void SetSigma(const double s) { m_sigma = s; }

    itk::RealTimeClock::TimeStampType GetTotalTime() const { return m_totalTime; }
    itk::RealTimeClock::TimeStampType GetMeanTime() const { return m_meanTime; }
    itk::SizeValueType GetEvaluations() const { return m_evaluations; }

    void GenerateOutputInformation() ITK_OVERRIDE {
        //std::cout <<  __PRETTY_FUNCTION__ << endl;
        Superclass::GenerateOutputInformation();
        const auto op = this->GetOutput();
        op->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
        op->SetNumberOfComponentsPerPixel(m_sequence->size());
    }

protected:
    SignalsFilter() {}
    ~SignalsFilter(){}

    void ThreadedGenerateData(const OutputImageRegionType &region, itk::ThreadIdType threadId) ITK_OVERRIDE {
        vector<itk::ImageRegionConstIterator<TImage>> inIters(m_model->nParameters());
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
                VectorXd parameters = m_model->Default();
                for (size_t i = 0; i < inIters.size(); i++) {
                    if (this->GetInput(i))
                        parameters[i] = inIters[i].Get();
                }
                VectorXcd allData = m_sequence->signal(m_model, parameters);
                if (m_sigma != 0.0) {
                    VectorXcd noise(m_sequence->size());
                    // Simple Box Muller transform
                    ArrayXd U = (ArrayXd::Random(m_sequence->size()) * 0.5) + 0.5;
                    ArrayXd V = (ArrayXd::Random(m_sequence->size()) * 0.5) + 0.5;
                    noise.real() = (m_sigma / M_SQRT2) * (-2. * U.log()).sqrt() * cos(2. * M_PI * V);
                    noise.imag() = (m_sigma / M_SQRT2) * (-2. * V.log()).sqrt() * sin(2. * M_PI * U);
                    allData += noise;
                }
                VectorXcf floatData = allData.cast<complex<float>>();
                itk::VariableLengthVector<complex<float>> dataVector(floatData.data(), m_sequence->size());
                outputIter.Set(dataVector);
                if (threadId == 0) {
                    m_clock.Stop();
                }
            }
            if (mask)
                ++maskIter;
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
        //std::cout <<  __PRETTY_FUNCTION__ << endl;
        itk::DataObject::Pointer output;
        if (idx < 1) {
            auto img = TCVImage::New();
            img->SetNumberOfComponentsPerPixel(m_sequence->size());
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

const string usage {
"Usage is: qisignal [options]\n\
\n\
Calculates multi-component DESPOT signals (mainly for testing purposes).\n\
The program will !*suppress for input (unless --no-!*suppress specified)\n\
\n\
All times (TR) are in SECONDS. All angles are in degrees.\n\
\n\""
};

//Read in all required files and data from cin
void ParseInput(vector<shared_ptr<SequenceBase>> &cs, vector<string> &names, bool suppress);
void ParseInput(vector<shared_ptr<SequenceBase>> &cs, vector<string> &names, bool suppress) {
    string type, path;
    if (!suppress) cout << "Enter output filename: " << flush;
    while (Read(cin, path) && (path != "END") && (path != "")) {
        names.push_back(path);
        cs.push_back(QI::ReadSequence(cin, !suppress));
        if (!suppress) cout << "Enter next filename (END to finish input): " << flush;
    }
}
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::OptionList opts(usage);
    QI::Option<int>    num_threads(4,'T',"threads","Use N threads (default=4, 0=hardware limit)", opts);
    QI::Switch         complex('x',"complex","Output complex valued signals", opts);
    QI::Option<float>  clamp(0.0,'c',"clamp","Clamp output between 0 and value", opts);
    QI::EnumOption     modelOpt("123",'1','M',"model","Choose number of components (1/2/3)", opts);
    QI::Option<int>    seed(-1,'s',"seed","Seed noise RNG with specific value", opts);
    QI::Option<float>  noise(0.0,'N',"noise","Add complex noise with std=value", opts);
    QI::ImageOption<QI::VolumeF> reference('r',"ref","Output space defined by reference volume", opts);
    QI::ImageOption<QI::VolumeF> mask('m', "mask", "Mask input with specified file", opts);
    QI::Option<string> outPrefix("", 'o', "out","Add a prefix to output filenames", opts);
    QI::Switch         suppress('n',"no-prompt","Suppress input prompts", opts);
    QI::Switch         verbose('v',"verbose","Print more information", opts);
    QI::Help help(opts);
    std::deque<std::string> nonopts = opts.parse(argc, argv);
    if (nonopts.size() != 0) {
        std::cerr << opts << std::endl;
        std::cerr << "Extraneous input argument provided." << std::endl;
        return EXIT_FAILURE;
    }

    std::shared_ptr<Model> model = nullptr;
    switch (*modelOpt) {
        case '1': model = make_shared<SCD>(); break;
        case '2': model = make_shared<MCD2>(); break;
        case '3': model = make_shared<MCD3>(); break;
    }
    if (*verbose) cout << "Using " << model->Name() << " model." << endl;

    if (seed.set()) {
        std::srand(*seed);
    } else {
        std::srand((unsigned int) time(0));
    }

    /***************************************************************************
     * Read in parameter files
     **************************************************************************/
    SignalsFilter::Pointer calcSignal = SignalsFilter::New();
    calcSignal->SetModel(model);
    calcSignal->SetSigma(*noise);
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(*num_threads);
    calcSignal->SetMask(*mask);
    if (*verbose) {
        auto monitor = QI::GenericMonitor::New();
        calcSignal->AddObserver(itk::ProgressEvent(), monitor);
    }
    if (!*suppress) cout << "Loading parameters." << endl;
    for (size_t i = 0; i < model->nParameters(); i++) {
        if (!*suppress) cout << "Enter path to " << model->ParameterNames()[i] << " file (blank for default value): " << flush;
        string filename;
        getline(cin, filename);
        if (filename != "") {
            if (*verbose) cout << "Opening " << filename << endl;
            QI::VolumeF::Pointer param = QI::ReadImage(filename);
            if (reference.set()) {
                if (*verbose) cout << "Resampling to reference" << endl;
                typedef itk::ResampleImageFilter<QI::VolumeF, QI::VolumeF, double> TResampler;
                typedef itk::LinearInterpolateImageFunction<QI::VolumeF, double> TInterp;
                typename TInterp::Pointer interp = TInterp::New();
                interp->SetInputImage(param);
                typename TResampler::Pointer resamp = TResampler::New();
                resamp->SetInput(param);
                resamp->SetInterpolator(interp);
                resamp->SetDefaultPixelValue(0.);
                resamp->SetOutputParametersFromImage(*reference);
                resamp->Update();
                QI::VolumeF::Pointer rparam = resamp->GetOutput();
                rparam->DisconnectPipeline();   
                calcSignal->SetInput(i, rparam);
            } else { 
                calcSignal->SetInput(i, param);
            }
        } else {
            if (*verbose) cout << "Using default value: " << model->Default()[i] << endl;
        }
    }

    /***************************************************************************
     * Set up sequences
     **************************************************************************/
    vector<shared_ptr<SequenceBase>> sequences;
    vector<string> filenames;
    ParseInput(sequences, filenames, *suppress);
    
    typedef itk::ClampImageFilter<QI::SeriesF, QI::SeriesF> TClamp;
    TClamp::Pointer clamp_filter = TClamp::New();
    if (clamp.set()) {
        clamp_filter->SetBounds(0, *clamp);
    }
    for (size_t i = 0; i < sequences.size(); i++) {
        if (*verbose) cout << "Generating sequence: " << endl << *(sequences[i]);
        calcSignal->SetSequence(sequences[i]);
        calcSignal->Update();
        if (*verbose) cout << "Mean evaluation time: " << calcSignal->GetMeanTime() << " s ( " << calcSignal->GetEvaluations() << " voxels)" << endl;
        if (*verbose) cout << "Converting to timeseries" << endl;
        QI::VectorToSeriesXF::Pointer vecTo4D = QI::VectorToSeriesXF::New();
        vecTo4D->SetInput(calcSignal->GetOutput());
        if (*verbose) cout << "Saving to filename: " << filenames[i] << endl;
        if (*complex) {
            QI::WriteImage<SeriesXF>(vecTo4D->GetOutput(), *outPrefix + filenames[i]);
        } else {
            auto abs = itk::ComplexToModulusImageFilter<QI::SeriesXF, QI::SeriesF>::New();
            abs->SetInput(vecTo4D->GetOutput());
            if (clamp.set()) {
                clamp_filter->SetInput(abs->GetOutput());
                QI::WriteImage<SeriesF>(clamp_filter->GetOutput(),filenames[i]);
            } else {
                QI::WriteImage<SeriesF>(abs->GetOutput(), *outPrefix + filenames[i]);
            }
        }
    }
    if (*verbose) cout << "Finished all sequences." << endl;
    return EXIT_SUCCESS;
}

