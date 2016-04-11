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
#include <getopt.h>
#include <exception>
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
#include "itkUnaryFunctorImageFilter.h"

#include "Filters/VectorToImageFilter.h"

#include "QI/Models/Model.h"
#include "QI/Sequences/Sequences.h"
#include "QI/Util.h"
#include "QI/ThreadPool.h"

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

	virtual void GenerateOutputInformation() override {
		//std::cout <<  __PRETTY_FUNCTION__ << endl;
		Superclass::GenerateOutputInformation();
		const auto op = this->GetOutput();
		op->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
		op->SetNumberOfComponentsPerPixel(m_sequence->size());
	}

protected:
	SignalsFilter() {}
	~SignalsFilter(){}

    virtual void ThreadedGenerateData(const OutputImageRegionType &region, itk::ThreadIdType threadId) override {
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

class EnsureFinite {
public:
    EnsureFinite() {};
    ~EnsureFinite() {};
    bool operator!=( const EnsureFinite & ) const { return false; }
    bool operator==( const EnsureFinite &other ) const { return !(*this != other); }
    inline complex<float> operator()( const complex<float> &v ) const {
        if (isfinite(real((proj(v))))) {
            return v;
        } else {
            cout << "Voxel value was: " << v << endl;
            return complex<float>(-1.0,0.0);
        }
    }
};

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: qisignal [options]\n\
\n\
Calculates multi-component DESPOT signals (mainly for testing purposes).\n\
The program will prompt for input (unless --no-prompt specified)\n\
\n\
All times (TR) are in SECONDS. All angles are in degrees.\n\
\n\
Options:\n\
    --help, -h        : Print this message.\n\
    --verbose, -v     : Print extra information.\n\
    --mask, -m file   : Only calculate inside the mask.\n\
    --out, -o path    : Add a prefix to the output filenames\n\
    --no-prompt, -n   : Don't print prompts for input.\n\
    --noise, -N val   : Add complex noise with std=val.\n\
    --seed, -S val    : Specify seed for noise.\n\
    --1, --2, --3     : Use 1, 2 or 3 component sequences (default 3).\n\
    --ref, -r FILE    : Output images in space defined by reference file.\n\
    --finite, -f      : Replace inf/NaN in output with zeros.\n\
    --complex, -x     : Output complex-valued signal.\n\
    --threads, -T N   : Use N threads (default=4, 0=hardware limit)\n"
};
shared_ptr<Model> model = make_shared<SCD>();
bool verbose = false, prompt = true, outputComplex = false;
string outPrefix = "";
size_t num_threads = 4;
const struct option long_opts[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"no-prompt", no_argument, 0, 'n'},
	{"noise", required_argument, 0, 'N'},
    {"seed", required_argument, 0, 'S'},
	{"1", no_argument, 0, '1'},
	{"2", no_argument, 0, '2'},
	{"3", no_argument, 0, '3'},
    {"ref", required_argument, 0, 'r'},
    {"finite", no_argument, 0, 'f'},
	{"complex", no_argument, 0, 'x'},
	{"threads", required_argument, 0, 'T'},
	{0, 0, 0, 0}
};
static const char *short_opts = "hvnN:S:m:o:123r:fxT:";
//******************************************************************************
#pragma mark Read in all required files and data from cin
//******************************************************************************
void ParseInput(vector<shared_ptr<SequenceBase>> &cs, vector<string> &names);
void ParseInput(vector<shared_ptr<SequenceBase>> &cs, vector<string> &names) {
    string type, path;
    if (prompt) cout << "Enter output filename: " << flush;
    while (Read(cin, path) && (path != "END") && (path != "")) {
        names.push_back(path);
        if (prompt) cout << "Enter sequence type: " << flush;
        Read(cin, type);
        if (type == "SPGR") {
			cs.push_back(make_shared<SPGRSimple>(prompt));
		} else if (type == "SPGR_ECHO") {
			cs.push_back(make_shared<SPGREcho>(prompt));
		} else if (type == "SPGR_FINITE") {
			cs.push_back(make_shared<SPGRFinite>(prompt));
		} else if (type == "SSFP") {
			cs.push_back(make_shared<SSFPSimple>(prompt));
        } else if (type == "SSFP_ECHO") {
            cs.push_back(make_shared<SSFPEcho>(prompt));
        } else if (type == "SSFP_ECHO_FLEX") {
            cs.push_back(make_shared<SSFPEchoFlex>(prompt));
        } else if (type == "SSFP_FINITE") {
			cs.push_back(make_shared<SSFPFinite>(prompt));
		} else if (type == "SSFP_GS") {
			cs.push_back(make_shared<SSFP_GS>(prompt));
		} else if (type == "IRSPGR") {
			cs.push_back(make_shared<IRSPGR>(prompt));
		} else if (type == "MPRAGE") {
			cs.push_back(make_shared<MPRAGE>(prompt));
		} else if (type == "AFI") {
			cs.push_back(make_shared<AFI>(prompt));
		} else if (type == "SPINECHO") {
			cs.push_back(make_shared<MultiEcho>(prompt));
		} else {
            QI_EXCEPTION("Unknown sequence type: " << type);
		}
        if (prompt) cout << "Enter next filename (END to finish input): " << flush;
	}
}
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	Eigen::initParallel();
	QI::ImageF::Pointer mask = ITK_NULLPTR, reference = ITK_NULLPTR;
    int indexptr = 0, c;
    double sigma = 0.;
    int seed = -1;
    bool ensure_finite = false;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
            case 'N': sigma = stod(optarg); break;
            case 'S': seed = stoi(optarg); break;
            case 'm':
                cout << "Reading mask file " << optarg << endl;
                mask = QI::ReadImage(optarg);
                break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case '1': model = make_shared<SCD>(); break;
			case '2': model = make_shared<MCD2>(); break;
			case '3': model = make_shared<MCD3>(); break;
            case 'r':
                if (verbose) cout << "Reading reference file: " << optarg << endl;
                reference = QI::ReadImage(optarg);
                break;
            case 'f':
                if (verbose) cout << "NaN/Inf will be set to 0 in output" << endl;
                ensure_finite = true;
                break;
			case 'x': outputComplex = true; break;
			case 'T':
                num_threads = stoi(optarg);
                if (num_threads == 0)
                    num_threads = std::thread::hardware_concurrency();
                break;
                if (verbose) cout << "Using " << num_threads << " threads" << endl;
			case 'h':
				cout << QI::GetVersion() << endl << usage << endl;
				return EXIT_SUCCESS;
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
			default:
				cout << "Unhandled option " << string(1, c) << endl;
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 0) {
		cerr << usage << endl << "Incorrect number of arguments." << endl;
		return EXIT_FAILURE;
	}  else if (prompt) {
		cout << "Starting qisignal" << endl;
		cout << "Run with -h switch to see usage" << endl;
	}
	if (verbose) cout << "Using " << model->Name() << " model." << endl;

    if (seed == -1) {
        std::srand((unsigned int) time(0));
    } else {
        std::srand(seed);
    }

	/***************************************************************************
	 * Read in parameter files
	 **************************************************************************/
	SignalsFilter::Pointer calcSignal = SignalsFilter::New();
	calcSignal->SetModel(model);
	calcSignal->SetSigma(sigma);
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(num_threads);
    calcSignal->SetMask(mask);
    if (verbose) {
        auto monitor = QI::GenericMonitor::New();
        calcSignal->AddObserver(itk::ProgressEvent(), monitor);
    }
	if (prompt) cout << "Loading parameters." << endl;
	for (size_t i = 0; i < model->nParameters(); i++) {
        if (prompt) cout << "Enter path to " << model->ParameterNames()[i] << " file (blank for default value): " << flush;
		string filename;
		getline(cin, filename);
        if (filename != "") {
            if (verbose) cout << "Opening " << filename << endl;
            QI::ImageF::Pointer param = QI::ReadImage(filename);
            if (reference) {
                if (verbose) cout << "Resampling to reference" << endl;
                typedef itk::ResampleImageFilter<QI::ImageF, QI::ImageF, double> TResampler;
                typedef itk::LinearInterpolateImageFunction<QI::ImageF, double> TInterp;
                typename TInterp::Pointer interp = TInterp::New();
                interp->SetInputImage(param);
                typename TResampler::Pointer resamp = TResampler::New();
                resamp->SetInput(param);
                resamp->SetInterpolator(interp);
                resamp->SetDefaultPixelValue(0.);
                resamp->SetOutputParametersFromImage(reference);
                resamp->Update();
                QI::ImageF::Pointer rparam = resamp->GetOutput();
                rparam->DisconnectPipeline();   
                calcSignal->SetInput(i, rparam);
            } else { 
                calcSignal->SetInput(i, param);
            }
        } else {
            if (verbose) cout << "Using default value: " << model->Default()[i] << endl;
        }
	}

	/***************************************************************************
	 * Set up sequences
	 **************************************************************************/
	vector<shared_ptr<SequenceBase>> sequences;
	vector<string> filenames;
    ParseInput(sequences, filenames);
    
    typedef itk::UnaryFunctorImageFilter<QI::TimeseriesXF, QI::TimeseriesXF, EnsureFinite> TFinite;
    TFinite::Pointer finite_filter = TFinite::New();
    QI::TimeseriesXF::Pointer out_series = ITK_NULLPTR;
	for (size_t i = 0; i < sequences.size(); i++) {
        if (verbose) cout << "Generating sequence: " << endl << *(sequences[i]);
		calcSignal->SetSequence(sequences[i]);
        calcSignal->Update();
        if (verbose) cout << "Mean evaluation time: " << calcSignal->GetMeanTime() << " s ( " << calcSignal->GetEvaluations() << " voxels)" << endl;
        if (verbose) cout << "Converting to timeseries" << endl;
        
        QI::VectorToTimeseriesXF::Pointer vecTo4D = QI::VectorToTimeseriesXF::New();
        vecTo4D->SetInput(calcSignal->GetOutput());
        if (ensure_finite) {
            if (verbose) cout << "Removing NaN/Inf" << endl;
            finite_filter->SetInput(vecTo4D->GetOutput());
            finite_filter->Update();
            out_series = finite_filter->GetOutput();
        } else {
            out_series = vecTo4D->GetOutput();
        }
        if (verbose) cout << "Saving to filename: " << filenames[i] << endl;
		if (outputComplex) {
			auto writer = QI::TimeseriesWriterXF::New();
            writer->SetInput(out_series);
			writer->SetFileName(filenames[i]);
			writer->Update();
		} else {
			auto writer = QI::TimeseriesWriterF::New();
			auto abs = itk::ComplexToModulusImageFilter<QI::TimeseriesXF, QI::TimeseriesF>::New();
            abs->SetInput(out_series);
			writer->SetInput(abs->GetOutput());
			writer->SetFileName(filenames[i]);
			writer->Update();
		}
	}
	if (verbose) cout << "Finished all sequences." << endl;
	return EXIT_SUCCESS;
}

