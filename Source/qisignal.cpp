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

#include "Filters/VectorToImageFilter.h"

#include "Model.h"
#include "Sequence.h"
#include "Util.h"
#include "ThreadPool.h"

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
    size_t m_nThreads = 1;
    itk::RealTimeClock::TimeStampType m_totalTime = 0.0;
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
    void SetPoolsize(const size_t n) { m_nThreads = n; }
    
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
    itk::RealTimeClock::TimeStampType GetMeanTime() const { return m_totalTime / m_evaluations; }
    itk::SizeValueType GetEvaluations() const { return m_evaluations; }

	virtual void GenerateOutputInformation() override {
		//std::cout <<  __PRETTY_FUNCTION__ << endl;
		Superclass::GenerateOutputInformation();
		const auto op = this->GetOutput();
		op->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
		op->SetNumberOfComponentsPerPixel(m_sequence->size());
		op->Allocate();
	}

	virtual void Update() override {
		//std::cout << __PRETTY_FUNCTION__ << endl;
		Superclass::Update();
	}

protected:
	SignalsFilter() {}
	~SignalsFilter(){}

    virtual void GenerateData() override {
		//std::cout <<  __PRETTY_FUNCTION__ << endl;
        TRegion region = this->GetInput(0)->GetLargestPossibleRegion();
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
        QI::ThreadPool threadPool(m_nThreads);
        itk::TimeProbe clock;
        clock.Start();
        m_evaluations = 0;
		while(!inIters[0].IsAtEnd()) {
			if (!mask || maskIter.Get()) {
                auto task = [=] {
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
                };
                threadPool.enqueue(task);
                ++m_evaluations;
			}
			if (mask)
				++maskIter;
			for (size_t i = 0; i < m_model->nParameters(); i++) {
                if (this->GetInput(i))
                    ++inIters[i];
			}
			++outputIter;
		}
        clock.Stop();
        m_totalTime = clock.GetTotal();
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


//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: mcsignal [options]\n\
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
	--complex, -x     : Output complex-valued signal.\n\
	--sequences, -M s : Use simple sequences (default).\n\
	            f     : Use Finite Pulse Length correction.\n\
	--threads, -T N   : Use N threads (default=1, 0=hardware limit)\n"
};
shared_ptr<Model> model = make_shared<SCD>();
bool verbose = false, prompt = true, finitesequences = false, outputComplex = false;
string outPrefix = "";
size_t num_threads = 1;
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
	{"complex", no_argument, 0, 'x'},
	{"sequences", no_argument, 0, 'M'},
	{"threads", required_argument, 0, 'T'},
	{0, 0, 0, 0}
};
static const char *short_opts = "hvnN:S:m:o:123xM:T:";
//******************************************************************************
#pragma mark Read in all required files and data from cin
//******************************************************************************
void parseInput(vector<shared_ptr<SequenceBase>> &cs, vector<string> &names);
void parseInput(vector<shared_ptr<SequenceBase>> &cs, vector<string> &names) {
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
		} else if (type == "SSFP_FINITE") {
			cs.push_back(make_shared<SSFPFinite>(prompt));
		} else if (type == "SSFP_GS") {
			cs.push_back(make_shared<SSFPEllipse>(prompt));
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
	QI::ImageF::Pointer mask = ITK_NULLPTR;
    int indexptr = 0, c;
    double sigma = 0.;
    int seed = -1;
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
			case 'x': outputComplex = true; break;
			case 'M':
				switch (*optarg) {
					case 's': finitesequences = false; if (prompt) cout << "Simple sequences selected." << endl; break;
					case 'f': finitesequences = true; if (prompt) cout << "Finite pulse correction selected." << endl; break;
					default:
						cout << "Unknown sequences type " << *optarg << endl;
						return EXIT_FAILURE;
						break;
				}
				break;
			case 'T':
                num_threads = stoi(optarg);
                if (num_threads == 0)
                    num_threads = std::thread::hardware_concurrency();
                break;
			case 'h':
				cout << usage << endl;
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
    calcSignal->SetPoolsize(num_threads);
    calcSignal->SetMask(mask);
	vector<QI::ImageReaderF::Pointer> pFiles(model->nParameters());
	if (prompt) cout << "Loading parameters." << endl;
	for (size_t i = 0; i < model->nParameters(); i++) {
        if (prompt) cout << "Enter path to " << model->ParameterNames()[i] << " file (blank for default value): " << flush;
		string filename;
		getline(cin, filename);
        if (filename != "") {
            if (verbose) cout << "Opening " << filename << endl;
            pFiles[i] = QI::ImageReaderF::New();
            pFiles[i]->SetFileName(filename);
            calcSignal->SetInput(i, pFiles[i]->GetOutput());
        } else {
            if (verbose) cout << "Using default value: " << model->Default()[i] << endl;
        }
	}

	/***************************************************************************
	 * Set up sequences
	 **************************************************************************/
	vector<shared_ptr<SequenceBase>> sequences;
	vector<string> filenames;
	parseInput(sequences, filenames);
	for (size_t i = 0; i < sequences.size(); i++) {
        if (verbose) cout << "Generating sequence: " << endl << *(sequences[i]);
		calcSignal->SetSequence(sequences[i]);
        calcSignal->Update();
        if (verbose) cout << "Mean evaluation time: " << calcSignal->GetMeanTime() << " s ( " << calcSignal->GetEvaluations() << " voxels)" << endl;
        if (verbose) cout << "Saving to filename: " << filenames[i] << endl;
		auto VecTo4D = QI::VectorToTimeseriesXF::New();
		VecTo4D->SetInput(calcSignal->GetOutput());
		if (outputComplex) {
			auto writer = QI::TimeseriesWriterXF::New();
			writer->SetInput(VecTo4D->GetOutput());
			writer->SetFileName(filenames[i]);
			writer->Update();
		} else {
			auto writer = QI::TimeseriesWriterF::New();
			auto abs = itk::ComplexToModulusImageFilter<QI::TimeseriesXF, QI::TimeseriesF>::New();
			abs->SetInput(VecTo4D->GetOutput());
			writer->SetInput(abs->GetOutput());
			writer->SetFileName(filenames[i]);
			writer->Update();
		}
	}
	if (verbose) cout << "Finished all sequences." << endl;
	return EXIT_SUCCESS;
}

