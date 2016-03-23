/*
 *  qmcdespot.cpp
 *
 *  Created by Tobias Wood on 03/06/2015.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <getopt.h>
#include <time.h>
#include <fstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "itkParticleSwarmOptimizer.h"
#include "itkLBFGSBOptimizer.h"
#include "itkTimeProbe.h"

#include "QI/Util.h"
#include "QI/Models/Model.h"
#include "QI/Sequences/Sequence.h"
#include "QI/RegionContraction.h"
#include "Filters/ApplyAlgorithmFilter.h"
#include "Filters/ReorderVectorFilter.h"

using namespace std;
using namespace Eigen;

/*
 * Read in all required files and data from cin
 */
void parseInput(shared_ptr<QI::SequenceGroup> seq, vector<typename QI::VectorImageF::Pointer> &images,
                bool flip, bool verbose, bool prompt);
void parseInput(shared_ptr<QI::SequenceGroup> seq, vector<typename QI::VectorImageF::Pointer> &images,
                bool flip, bool verbose, bool prompt)
{
    string type, path;
    if (prompt) cout << "Enter input filename: " << flush;
    while (QI::Read(cin, path) && (path != "END") && (path != "")) {
        auto file = QI::TimeseriesReaderF::New();
        auto data = QI::TimeseriesToVectorF::New();
        QI::VectorImageF::Pointer image;
        if (verbose) cout << "Reading file: " << path << endl;
        file->SetFileName(path);
        data->SetInput(file->GetOutput());
        data->Update();
        image = data->GetOutput();
        if (prompt) cout << "Enter sequence type (SPGR/SSFP): " << flush;
        QI::Read(cin, type);
        if (type == "SPGR") {
            seq->addSequence(make_shared<QI::SPGRSimple>(prompt));
        } else if (type == "SPGR_ECHO") {
            seq->addSequence(make_shared<QI::SPGREcho>(prompt));
        } else if (type == "SPGR_FINITE") {
            seq->addSequence(make_shared<QI::SPGRFinite>(prompt));
        } else if (type == "SSFP") {
            seq->addSequence(make_shared<QI::SSFPSimple>(prompt));
        } else if (type == "SSFP_ECHO") {
            seq->addSequence(make_shared<QI::SSFPEcho>(prompt));
        } else if (type == "SSFP_ECHO_FLEX") {
            seq->addSequence(make_shared<QI::SSFPEchoFlex>(prompt));
        } else if (type == "SSFP_FINITE") {
            seq->addSequence(make_shared<QI::SSFPFinite>(prompt));
        } else {
            QI_EXCEPTION("Unknown sequence type: " << type);
        }
        image->DisconnectPipeline(); // This step is really important.
        images.push_back(image);
        if (prompt) cout << "Enter next filename (END to finish input): " << flush;
    }
}

namespace itk {
class MCDCostFunction : public SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef MCDCostFunction           Self;
    typedef SingleValuedCostFunction  Superclass;
    typedef SmartPointer<Self>        Pointer;
    typedef SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(CostFunction, SingleValuedCostfunction);

    shared_ptr<QI::SequenceGroup> m_sequence;
    ArrayXd m_data, m_weights;
    shared_ptr<QI::Model> m_model;

    unsigned int GetNumberOfParameters(void) const override { return m_model->nParameters(); } // itk::CostFunction

    ArrayXd residuals(const ParametersType &p1) const {
        int N = m_model->nParameters();
        ArrayXd p2(N);
        for (int i = 0; i < N; i++)
            p2[i] = p1[i];
        ArrayXcd s = m_sequence->signal(m_model, p2);
        return (s.abs() - m_data);
    }

    MeasureType GetValue(const ParametersType &p1) const override {
        return (residuals(p1) * m_weights).square().sum();
    }

    void GetDerivative(const ParametersType &p, DerivativeType &df) const override {
        using std::sqrt;
        using std::abs;
        double eps = sqrt(numeric_limits<double>::epsilon());
        df.SetSize(p.Size());
        ParametersType x = p;
        for (int i = 0; i < p.Size(); ++i) {
            double h = eps * abs(x[i]);
            if (h == 0.) {
                h = eps;
            }
            double _x = x[i];
            x[i] += h;
            double v2 = GetValue(x);
            x[i] -= 2*h;
            double v1 = GetValue(x);
            df[i] = (v2 - v1)/(2*h);
            x[i] = _x;
        }
        //cout << "p " << p << endl;
        //cout << "x " << x.transpose() << endl;
        //cout << "e " << (eps * x).transpose() << endl;
        //cout << "df " << df << endl;
        //exit(0);
    }

    protected:
    MCDCostFunction(){};
    ~MCDCostFunction(){};

    private:
    MCDCostFunction(const Self &); //purposely not implemented
    void operator = (const Self &); //purposely not implemented
};
}

class MCDAlgo : public Algorithm<double> {
protected:
    ArrayXXd m_bounds;
    shared_ptr<QI::Model> m_model = nullptr;
    shared_ptr<QI::SequenceGroup> m_sequence = nullptr;
    QI::FieldStrength m_tesla = QI::FieldStrength::Three;
    int m_iterations = 0;

public:
    MCDAlgo(shared_ptr<QI::Model>&m, ArrayXXd &b,
            shared_ptr<QI::SequenceGroup> s, int mi) :
        m_model(m), m_bounds(b), m_sequence(s), m_iterations(mi)
    {}

    size_t numInputs() const override  { return m_sequence->count(); }
    size_t numOutputs() const override { return m_model->nParameters(); }
    size_t dataSize() const override   { return m_sequence->size(); }

    void setModel(shared_ptr<QI::Model> &m) { m_model = m; }
    void setSequence(shared_ptr<QI::SequenceGroup> &s) { m_sequence = s; }
    void setBounds(ArrayXXd &b) { m_bounds = b; }
    void setIterations(const int i) { m_iterations = i; }
};

class BFGSAlgo : public MCDAlgo {
    using MCDAlgo::MCDAlgo;

private:
    ArrayXd m_start;

public:
    void setStart(ArrayXd &s) { m_start = s; }

    size_t numConsts() const override  { return 2; }
    virtual TArray defaultConsts() override {
        // f0, B1
        TArray def(4);
        def << NAN, 1.;
        return def;
    }

    virtual void apply(const TInput &data, const TArray &inputs,
                       TArray &outputs, TArray &resids, TIterations &its) const override
    {
        typedef itk::LBFGSBOptimizer TOpt;
        TOpt::Pointer optimizer = TOpt::New();
        auto cost = itk::MCDCostFunction::New();
        double f0 = inputs[0];
        double B1 = inputs[1];
        ArrayXXd localBounds = m_bounds;
        if (isfinite(f0)) { // We have an f0 map, add it to the fitting bounds
            localBounds.row(m_model->ParameterIndex("f0")) += f0;
            cost->m_weights = m_sequence->weights(f0);
        } else {
            cost->m_weights = ArrayXd::Ones(m_sequence->size());
        }
        localBounds.row(m_model->ParameterIndex("B1")).setConstant(B1);
        int N = m_model->nParameters();
        TOpt::ParametersType start(N);
        TOpt::BoundSelectionType select(N);
        TOpt::BoundValueType lower(N), upper(N);
        for (int i = 0; i < N; i++) {
            lower[i] = localBounds(i, 0);
            upper[i] = localBounds(i, 1);
            select[i] = 2; // Upper and lower bounds
            start[i] = m_start[i];
        }
        cost->m_data = data;
        cost->m_sequence = m_sequence;
        cost->m_model = m_model;
        //optimizer->DebugOn();
        optimizer->SetCostFunction(cost);
        optimizer->SetBoundSelection(select);
        optimizer->SetLowerBound(lower);
        optimizer->SetUpperBound(upper);
        optimizer->SetCostFunctionConvergenceFactor(1.e3);
        optimizer->SetProjectedGradientTolerance(1e-10);
        optimizer->SetMaximumNumberOfIterations(m_iterations);
        optimizer->SetMaximumNumberOfEvaluations(999999);
        optimizer->SetMaximumNumberOfCorrections(20);
        optimizer->SetInitialPosition(start);
        //optimizer->SetTrace(true);
        optimizer->StartOptimization();
        TOpt::ParametersType final = optimizer->GetCurrentPosition();
        //outputs[0] = m_PDscale;
        for (int i = 0; i < m_model->nParameters(); i++) {
            outputs[i] = final[i];
        }
        //outputs[m_model->ParameterIndex("f0")] = f0;
        //outputs[m_model->ParameterIndex("B1")] = B1;
        resids = cost->residuals(final);
        its = optimizer->GetCurrentIteration();
    }
};

class MCDSRCFunctor {
    public:
        const shared_ptr<QI::SequenceGroup> m_sequence;
        const ArrayXd m_data, m_weights;
        const shared_ptr<QI::Model> m_model;

        MCDSRCFunctor(shared_ptr<QI::Model> m,shared_ptr<QI::SequenceGroup> s, const ArrayXd &d, const ArrayXd &w) :
            m_sequence(s), m_data(d), m_model(m), m_weights(w)
        {
            assert(static_cast<size_t>(m_data.rows()) == m_sequence->size());
        }

        int inputs() const { return m_model->nParameters(); }
        int values() const { return m_sequence->size(); }

        const bool constraint(const VectorXd &params) const {
            return m_model->ValidParameters(params);
        }

        ArrayXd residuals(const Ref<VectorXd> &params) const {
            const ArrayXd s = (m_sequence->signal(m_model, params)).abs();
            return m_data - s;
        }
        double operator()(const Ref<VectorXd> &params) const {
            eigen_assert(diffs.size() == values());
            return (residuals(params) * m_weights).square().sum();
        }
};

class SRCAlgo : public MCDAlgo {
    using MCDAlgo::MCDAlgo;

    private:
        size_t m_samples = 5000, m_retain = 50;
        bool m_gauss = true;

    public:
        void setGauss(bool g) { m_gauss = g; }

        size_t numConsts() const override  { return 2; }
        virtual TArray defaultConsts() override {
			// f0, B1
			TArray def = TArray::Ones(2);
            def[0] = NAN;
			return def;
		}

		virtual void apply(const TInput &data, const TArray &inputs,
                           TArray &outputs, TArray &resids, TIterations &its) const override
		{
			ArrayXd thresh(m_model->nParameters()); thresh.setConstant(0.05);
			double f0 = inputs[0];
			double B1 = inputs[1];
            ArrayXXd localBounds = m_bounds;
            ArrayXd weights = ArrayXd::Ones(m_sequence->size());
            if (isfinite(f0)) { // We have an f0 map, add it to the fitting bounds
                localBounds.row(m_model->ParameterIndex("f0")) += f0;
                weights = m_sequence->weights(f0);
            }
            localBounds.row(m_model->ParameterIndex("B1")).setConstant(B1);
            MCDSRCFunctor func(m_model, m_sequence, data, weights);
            QI::RegionContraction<MCDSRCFunctor> rc(func, localBounds, thresh, m_samples, m_retain, m_iterations, 0.02, m_gauss, false);
            rc.optimise(outputs);
            //outputs(m_model->nParameters() - 1) = rc.contractions();
            //outputs(0) = static_cast<int>(rc.status());
            resids = func.residuals(outputs);
            its = rc.contractions();
		}
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();

	auto tesla = QI::FieldStrength::Three;
	int start_slice = 0, stop_slice = 0;
    int max_its = 4, num_threads = 4;
	int verbose = false, prompt = true, all_residuals = false, flipData = false;
    string outPrefix;
    enum class Algos { SRC, GRC, BFGS };
    Algos which_algo = Algos::GRC;

	QI::ImageReaderF::Pointer mask, B1, f0 = ITK_NULLPTR;
	shared_ptr<QI::Model> model = make_shared<QI::MCD3>();
	typedef itk::VectorImage<float, 2> VectorSliceF;
    typedef itk::ApplyAlgorithmFilter<MCDAlgo> TMCDFilter;
	auto apply = TMCDFilter::New();
    const string usage {
"Usage is: qimcdespot [options]\n\
\n\
The program will prompt for input (unless --no-prompt specified)\n\
\n\
All times (TR) are in SECONDS. All angles are in degrees.\n\
\n\
Options:\n\
    --help, -h        : Print this message\n\
    --verbose, -v     : Print more information\n\
    --no-prompt, -n   : Don't print prompts for input\n\
    --mask, -m file   : Mask input with specified file\n\
    --out, -o path    : Add a prefix to the output filenames\n\
    --model, -M 1     : Use 1 component model\n\
                2     : Use 2 component model\n\
                2nex  : Use 2 component, no exchange model\n\
                3     : Use 3 component model (default)\n\
                3_f0  : Use 3 components with a myelin resonance frequency\n\
                3nex  : Use 3 component, no exchange model\n\
    --f0, -f file     : Use f0 Map file (in Hertz)\n\
    --B1, -b file     : B1 Map file (ratio)\n\
    --start, -s n     : Only start processing at slice n.\n\
    --stop, -p n      : Finish at slice n-1\n\
    --scale, -S       : Normalise signals to mean\n\
    --algo, -a S      : Use Uniform distribution for Region Contraction\n\
               G      : Use Gaussian distribution for RC (default)\n\
               b      : Use BFGS algorithm\n\
    --iters, -i N     : Specify maximum number of iterations (default 4)\n\
    --flip, -F        : Data order is phase then flip-angle (default opposite)\n\
    --tesla, -t 3     : Boundaries suitable for 3T (default)\n\
                7     : Boundaries suitable for 7T \n\
                u     : User specified boundaries from stdin\n\
    --resids, -r      : Write out per flip-angle residuals\n\
    --threads, -T N   : Use N threads (default=4, 0=hardware limit)\n"
	};

	const struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"verbose", no_argument, 0, 'v'},
		{"mask", required_argument, 0, 'm'},
		{"out", required_argument, 0, 'o'},
		{"f0", required_argument, 0, 'f'},
		{"B1", required_argument, 0, 'b'},
		{"start", required_argument, 0, 's'},
		{"stop", required_argument, 0, 'p'},
        {"scale", no_argument, 0, 'S'},
        {"algo", required_argument, 0, 'a'},
        {"iterations", required_argument, 0, 'i'},
        {"flip", no_argument, 0, 'F'},
		{"tesla", required_argument, 0, 't'},
		{"resids", no_argument, 0, 'r'},
		{"threads", required_argument, 0, 'T'},
		{"no-prompt", no_argument, 0, 'n'},
        {"model", required_argument, 0, 'M'},
		{0, 0, 0, 0}
	};
    const char* short_options = "hvm:o:f:b:s:p:Sa:t:FT:rnM:i:j:";

	// Deal with these options in first pass to ensure the correct model is selected
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
            case 'M': {
                string choose_model(optarg);
                if (choose_model == "1") { model = make_shared<QI::SCD>(); }
                else if (choose_model == "2") { model = make_shared<QI::MCD2>(); }
                else if (choose_model == "2nex") { model = make_shared<QI::MCD2_NoEx>(); }
                else if (choose_model == "3") { model = make_shared<QI::MCD3>(); }
                else if (choose_model == "3_f0") { model = make_shared<QI::MCD3_f0>(); }
                else if (choose_model == "3nex") { model = make_shared<QI::MCD3_NoEx>(); }
            } break;
			default:
				break;
		}
	}
	// Now reset and do a second pass
	optind = 1;
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
            case 'v': case 'n': case 'M': break; // Already handled
			case 'm':
				if (verbose) cout << "Reading mask file " << optarg << endl;
				mask = QI::ImageReaderF::New();
				mask->SetFileName(optarg);
				break;
			case 'o':
				outPrefix = optarg;
				if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'f':
				if (verbose) cout << "Reading f0 file: " << optarg << endl;
				f0 = QI::ImageReaderF::New();
				f0->SetFileName(optarg);
				break;
			case 'b':
				if (verbose) cout << "Reading B1 file: " << optarg << endl;
				B1 = QI::ImageReaderF::New();
				B1->SetFileName(optarg);
				break;
			case 's': start_slice = atoi(optarg); break;
			case 'p': stop_slice = atoi(optarg); break;
            case 'S':
                if (verbose) cout << "Mean scaling selected." << endl;
                model->setScaleToMean(true);
                apply->SetScaleToMean(true);
                break;
            case 'a':
                switch (*optarg) {
                case 'S': which_algo = Algos::SRC; break;
                case 'G': which_algo = Algos::GRC; break;
                case 'b': which_algo = Algos::BFGS; break;
                default:
                    cerr << "Unknown algorithm type " << *optarg << endl;
                    return EXIT_FAILURE;
                    break;
                } break;
            case 'F': flipData = true; if (verbose) cout << "Data order is phase, then flip-angle" << endl; break;
			case 'T': 
                num_threads = stoi(optarg);
                if (num_threads == 0)
                    num_threads = std::thread::hardware_concurrency();
                break;
			case 't':
				switch (*optarg) {
					case '3': tesla = QI::FieldStrength::Three; break;
					case '7': tesla = QI::FieldStrength::Seven; break;
					case 'u': tesla = QI::FieldStrength::User; break;
					default:
						cerr << "Unknown boundaries type " << *optarg << endl;
						return EXIT_FAILURE;
						break;
                } break;
            case 'i': max_its = atoi(optarg); break;
			case 'r': all_residuals = true; break;
			case 'h':
				cout << QI::GetVersion() << endl << usage << endl;
				return EXIT_SUCCESS;
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
			case 0: break; // Just a flag
			default:
				cout << "Unhandled option " << string(1, c) << endl;
				return EXIT_FAILURE;
		}
	}
    if ((argc - optind) != 0) {
		cerr << usage << endl << "Incorrect number of arguments." << endl;
		return EXIT_FAILURE;
	} else if (prompt) {
		cout << "Starting qimcdespot" << endl;
		cout << "Run with -h switch to see usage" << endl;
	}
    Array2d f0Bandwidth;

	shared_ptr<QI::SequenceGroup> sequences = make_shared<QI::SequenceGroup>();
	// Build a Functor here so we can query number of parameters etc.
	if (verbose) cout << "Using " << model->Name() << " model." << endl;
    vector<QI::VectorImageF::Pointer> images;
    parseInput(sequences, images, flipData, verbose, prompt);

    ArrayXXd bounds = model->Bounds(tesla);
    ArrayXd start = model->Default(tesla);
    if (tesla == QI::FieldStrength::User) {
        ArrayXd temp;
        if (prompt) cout << "Enter lower bounds" << endl;
        QI::ReadArray(cin, temp);
        bounds.col(0) = temp;
        if (prompt) cout << "Enter upper bounds" << endl;
        QI::ReadArray(cin, temp);
        bounds.col(1) = temp;
        if (which_algo == Algos::BFGS) {
            if (prompt) cout << "Enter start point" << endl;
            QI::ReadArray(cin, start);
        }
    }
    switch (which_algo) {
    case Algos::SRC: {
        if (verbose) cout << "Using SRC algorithm" << endl;
        shared_ptr<SRCAlgo> algo = make_shared<SRCAlgo>(model, bounds, sequences, max_its);
        algo->setGauss(false);
        apply->SetAlgorithm(algo);
    } break;
    case Algos::GRC: {
        if (verbose) cout << "Using GRC algorithm" << endl;
        shared_ptr<SRCAlgo> algo = make_shared<SRCAlgo>(model, bounds, sequences, max_its);
        algo->setGauss(true);
        apply->SetAlgorithm(algo);
    } break;
    case Algos::BFGS: {
        if (verbose) cout << "Using BFGS algorithm" << endl;
        itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);
        shared_ptr<BFGSAlgo> algo = make_shared<BFGSAlgo>(model, bounds, sequences, max_its);
        algo->setStart(start);
        apply->SetAlgorithm(algo);
    } break;
    }
    apply->SetVerbose(verbose);
	apply->SetSlices(start_slice, stop_slice);
    apply->SetPoolsize(num_threads);
    for (int i = 0; i < images.size(); i++) {
        apply->SetInput(i, images[i]);
    }
	if (f0) {
		f0->Update();
		apply->SetConst(0, f0->GetOutput());
	}
	if (B1) {
		B1->Update();
		apply->SetConst(1, B1->GetOutput());
	}
	if (mask) {
		mask->Update();
		apply->SetMask(mask->GetOutput());
	}

    // Need this here so the bounds.txt file will have the correct prefix
    outPrefix = outPrefix + model->Name() + "_";
    if (verbose) {
        cout << *sequences;
        cout << "Bounds:" << endl <<  bounds.transpose() << endl;
        if (which_algo == Algos::BFGS)
            cout << "Start: " << endl << start.transpose() << endl;
        ofstream boundsFile(outPrefix + "bounds.txt");
        boundsFile << "Names: ";
        for (size_t p = 0; p < model->nParameters(); p++) {
            boundsFile << model->ParameterNames()[p] << "\t";
        }
        boundsFile << endl << "Bounds: " << endl << bounds.transpose() << endl;
        if (which_algo == Algos::BFGS)
            boundsFile << "Start: " << endl << start.transpose() << endl;
        boundsFile.close();
    }

    if (verbose) {
        cout << "Processing" << endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        cout << "Elapsed time was " << apply->GetTotalTime() << "s" << endl;
        cout << "Mean time per voxel was " << apply->GetMeanTime() << "s " << endl;
        cout << "Writing results files." << endl;
    }
	for (int i = 0; i < model->nParameters(); i++) {
        QI::WriteImage(apply->GetOutput(i), outPrefix + model->ParameterNames()[i] + QI::OutExt());
	}
	QI::WriteResiduals(apply->GetResidOutput(), outPrefix, all_residuals, apply->GetOutput(0));
    QI::WriteImage(apply->GetIterationsOutput(), outPrefix + "iterations" + QI::OutExt());
	return EXIT_SUCCESS;
}

