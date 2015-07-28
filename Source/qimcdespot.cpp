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

#include "Util.h"
#include "Filters/ApplyAlgorithmSliceBySliceFilter.h"
#include "Filters/ReorderVectorFilter.h"
#include "Model.h"
#include "Sequence.h"
#include "RegionContraction.h"

using namespace std;
using namespace Eigen;

/*
 * Read in all required files and data from cin
 */
void parseInput(shared_ptr<SequenceGroup> seq,
                vector<typename QI::ReadTimeseriesF::Pointer> &files,
                vector<typename QI::TimeseriesToVectorF::Pointer> &data,
                vector<typename QI::ReorderF::Pointer> &order,
                Array2d &f0Bandwidth, bool finite, bool flip, bool verbose, bool prompt);
void parseInput(shared_ptr<SequenceGroup> seq,
                vector<typename QI::ReadTimeseriesF::Pointer> &files,
                vector<typename QI::TimeseriesToVectorF::Pointer> &data,
                vector<typename QI::ReorderF::Pointer> &order,
                Array2d &f0Bandwidth, bool finite, bool flip, bool verbose, bool prompt)
{
	string type, path;
	if (verbose && finite) cout << "Using finite pulse-width sequences." << endl;
	if (prompt) cout << "Specify next image type (SPGR/SSFP): " << flush;
	f0Bandwidth = Array2d::Zero();
	while (QI::Read(cin, type) && (type != "END") && (type != "")) {
		if (type != "SPGR" && type != "SSFP") {
			throw(std::runtime_error("Unknown signal type: " + type));
		}
		if (prompt) cout << "Enter image path: " << flush;
		QI::Read(cin, path);
		files.push_back(QI::ReadTimeseriesF::New());
		files.back()->SetFileName(path);
		data.push_back(QI::TimeseriesToVectorF::New());
		data.back()->SetInput(files.back()->GetOutput());
		order.push_back(QI::ReorderF::New());
		order.back()->SetInput(data.back()->GetOutput());
		if (type == "SPGR") {
			if (finite) {
				seq->addSequence(make_shared<SPGRFinite>(prompt));
			} else {
				seq->addSequence(make_shared<SPGRSimple>(prompt));
			}
		} else if (type == "SSFP") {
			shared_ptr<SSFPSimple> s;
			if (finite) {
				s = make_shared<SSFPFinite>(prompt);
			} else {
				s = make_shared<SSFPSimple>(prompt);
			}
            f0Bandwidth(1) = 0.5 / s->TR();
            if (s->isSymmetric()) {
                f0Bandwidth(0) = 0.;
            } else {
                f0Bandwidth(0) = -f0Bandwidth(1);
            }
			seq->addSequence(s);
			if (flip)
				order.back()->SetStride(s->phases());
		}
		if (verbose) cout << "Reading file: " << path << endl;
		order.back()->Update();
		// Print message ready for next loop
		if (prompt) cout << "Specify next image type (SPGR/SSFP, END to finish input): " << flush;
	}
}

class MCDFunctor : public DenseFunctor<double> {
	public:
		const shared_ptr<SequenceBase> m_sequence;
		const ArrayXd m_data;
		const shared_ptr<Model> m_model;

        MCDFunctor(shared_ptr<Model> m,shared_ptr<SequenceBase> s, const ArrayXd &d) :
			DenseFunctor<double>(m->nParameters(), s->size()),
			m_sequence(s), m_data(d), m_model(m)
		{
			assert(static_cast<size_t>(m_data.rows()) == values());
		}

		const bool constraint(const VectorXd &params) const {
			return m_model->ValidParameters(params);
		}

		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
			eigen_assert(diffs.size() == values());
			ArrayXcd s = m_sequence->signal(m_model, params);
			diffs = m_data - s.abs();
			//cout << __PRETTY_FUNCTION__ << std::endl;
			//cout << "d " << m_data.transpose() << endl;
			//cout << "s " << s.abs().transpose() << endl;
			return 0;
		}
};

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

    shared_ptr<SequenceGroup> m_sequence;
    ArrayXd m_data;
    shared_ptr<Model> m_model;
    double m_PD, m_B1, m_f0;

    unsigned int GetNumberOfParameters(void) const { return m_model->nParameters() - 3; } // itk::CostFunction

    ArrayXd residuals(const ParametersType &p1) const {
        int N = m_model->nParameters();
        ArrayXd p2(N);
        p2[0] = m_PD;
        for (int i = 1; i < N - 2; i++)
            p2[i] = p1[i-1];
        p2[N-2] = m_f0;
        p2[N-1] = m_B1;
        ArrayXcd s = m_sequence->signal(m_model, p2);
        return (s.abs() - m_data);
    }

    MeasureType GetValue(const ParametersType &p1) const {
        return (residuals(p1) * m_sequence->weights(m_f0)).square().sum();
    }

    void GetDerivative(const ParametersType &p, DerivativeType &df) const {
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
    shared_ptr<Model> m_model = nullptr;
    shared_ptr<SequenceGroup> m_sequence = nullptr;
    FieldStrength m_tesla = FieldStrength::Three;
    double m_PDscale = 0.;

public:
    size_t numInputs() const override  { return m_sequence->count(); }
    size_t numOutputs() const override { return m_model->nParameters(); }
    size_t dataSize() const override   { return m_sequence->size(); }

    void setModel(shared_ptr<Model> &m) { m_model = m; }
    void setSequence(shared_ptr<SequenceGroup> &s) { m_sequence = s; }
    void setBounds(ArrayXXd &b) { m_bounds = b; }
    void setPDScale(double s) { m_PDscale = s; }
    void setTesla(FieldStrength t) { m_tesla = t; }
};

class LBFGSBAlgo : public MCDAlgo {
private:
    ArrayXd m_start;

public:
    void setStart(ArrayXd &s) { m_start = s; }

    size_t numConsts() const override  { return 5; }
    virtual TArray defaultConsts() {
        // f0, B1 T1, T2
        TArray def(4);
        def << 0., 1., 1., 1.;
        return def;
    }

    virtual void apply(const TInput &data, const TArray &inputs,
                       TArray &outputs, TArray &resids, TIterations &its) const override
    {
        double f0 = inputs[0];
        double B1 = inputs[1];
        double T1 = inputs[2];
        double T2 = inputs[3];
        ArrayXXd localBounds = m_bounds;
        ArrayXd localStart = m_model->Start(m_tesla, T1, T2);
        localBounds.row(0).setConstant(1.);
        localBounds.row(m_model->nParameters() - 2).setConstant(f0);
        localBounds.row(m_model->nParameters() - 1).setConstant(B1);

        typedef itk::LBFGSBOptimizer TOpt;
        int N = m_model->nParameters() - 3;
        TOpt::ParametersType start(N);
        TOpt::BoundSelectionType select(N);
        TOpt::BoundValueType lower(N), upper(N);
        for (int i = 0; i < N; i++) {
            lower[i] = localBounds(i+1, 0);
            upper[i] = localBounds(i+1, 1);
            select[i] = 2; // Upper and lower bounds
            start[i] = localStart[i+1];
        }
        //initP[0] = upper[0];
        //initP[2] = lower[2];
        TOpt::Pointer optimizer = TOpt::New();
        optimizer->SetCostFunctionConvergenceFactor(1.e3);
        optimizer->SetProjectedGradientTolerance(1e-10);
        optimizer->SetMaximumNumberOfIterations(150);
        optimizer->SetMaximumNumberOfEvaluations(999999);
        optimizer->SetMaximumNumberOfCorrections(20);
        auto cost = itk::MCDCostFunction::New();
        cost->m_data = data;
        cost->m_sequence = m_sequence;
        cost->m_model = m_model;
        cost->m_PD = 1.0;
        cost->m_f0 = f0;
        cost->m_B1 = B1;
        //optimizer->DebugOn();
        optimizer->SetCostFunction(cost);
        optimizer->SetBoundSelection(select);
        optimizer->SetLowerBound(lower);
        optimizer->SetUpperBound(upper);
        optimizer->DebugOff();
        optimizer->SetGlobalWarningDisplay(false);
        optimizer->SetInitialPosition(start);
        optimizer->StartOptimization();
        //cout << "Stop: " << optimizer->GetStopConditionDescription() << endl;
        TOpt::ParametersType final = optimizer->GetCurrentPosition();
        outputs[0] = m_PDscale;
        for (int i = 1; i < m_model->nParameters()-2; i++) {
            outputs[i] = final[i-1];
        }
        outputs[m_model->nParameters()-2] = f0;
        outputs[m_model->nParameters()-1] = B1;
        resids = cost->residuals(final);
        its = optimizer->GetCurrentIteration();
    }
};

class SRCAlgo : public MCDAlgo {
	private:
        size_t m_samples = 5000, m_retain = 50, m_contractions = 21;
        bool m_gauss = true;

	public:
		void setRCPars(size_t c, size_t s, size_t r) { m_contractions = c; m_samples = s; m_retain = r; }
		void setGauss(bool g) { m_gauss = g; }

		size_t numConsts() const override  { return 2; }
		virtual TArray defaultConsts() {
			// f0, B1
			TArray def = TArray::Ones(2);
            def[0] = NAN;
			return def;
		}

		virtual void apply(const TInput &data, const TArray &inputs,
                           TArray &outputs, TArray &resids, TIterations &its) const override
		{
			ArrayXd thresh(m_model->nParameters()); thresh.setConstant(0.05);
			ArrayXd weights = ArrayXd::Ones(m_sequence->size());
			double f0 = inputs[0];
			double B1 = inputs[1];
            ArrayXXd localBounds = m_bounds;
			if (!std::isnan(f0)) {
				localBounds.row(m_model->nParameters() - 2).setConstant(f0);
				weights = m_sequence->weights(f0);
			}
			if (m_PDscale == 0.) { // Then we need a reasonable fitting range for PD
				localBounds(0, 0) = 0.;
				localBounds(0, 1) = data.array().maxCoeff() * 25;
			} else {
				localBounds.row(0).setConstant(m_PDscale);
			}
			localBounds.row(m_model->nParameters() - 1).setConstant(B1);
            MCDFunctor func(m_model, m_sequence, data);
            RegionContraction<MCDFunctor> rc(func, localBounds, weights, thresh,
                                            m_samples, m_retain, m_contractions, 0.02, m_gauss, false);
            rc.optimise(outputs);
            //outputs(m_model->nParameters() - 1) = rc.contractions();
            //outputs(0) = static_cast<int>(rc.status());
            resids = rc.residuals();
            its = rc.contractions();
		}
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();

	auto tesla = FieldStrength::Three;
	int start_slice = 0, stop_slice = 0;
    int nargs = 0;
	int verbose = false, prompt = true, all_residuals = false,
        fitFinite = false, flipData = false;
	string outPrefix;
    double PDScale = 1.;
    enum class Algos { SRC, GRC, LBFGSB };
    Algos algo = Algos::SRC;

	QI::ReadImageF::Pointer mask, B1, f0 = ITK_NULLPTR;
	shared_ptr<Model> model = make_shared<MCD3>();
	typedef itk::VectorImage<float, 2> VectorSliceF;
    typedef itk::ApplyAlgorithmSliceBySliceFilter<MCDAlgo> TMCDFilter;
	auto applySlices = TMCDFilter::New();

	const string usage {
	"Usage is: mcdespot [options]\n\
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
                    3nex  : Use 3 component, no exchange model\n\
		--f0, -f file     : Use f0 Map file (in Hertz)\n\
		--B1, -b file     : B1 Map file (ratio)\n\
		--start, -s n     : Only start processing at slice n.\n\
		--stop, -p n      : Finish at slice n-1\n\
		--scale, -S MEAN  : Normalise signals to mean (default)\n\
					0     : Fit a scaling factor/proton density\n\
					val   : Fix PD to val\n\
        --algo, -a S      : Use Uniform distribution for Region Contraction\n\
                   G      : Use Gaussian distribution for RC (default)\n\
                   L      : Use LBFGS algorithm\n\
		--flip, -F        : Data order is phase, then flip-angle (default opposite)\n\
		--tesla, -t 3     : Boundaries suitable for 3T (default)\n\
					7     : Boundaries suitable for 7T \n\
					u     : User specified boundaries from stdin\n\
		--finite          : Use Finite Pulse Length correction\n\
		--contract, -c n  : Read contraction settings from stdin (Will prompt)\n\
		--resids, -r      : Write out per flip-angle residuals\n\
		--threads, -T N   : Use N threads (default=hardware limit)\n"
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
		{"scale", required_argument, 0, 'S'},
        {"algo", required_argument, 0, 'a'},
		{"flip", required_argument, 0, 'F'},
		{"tesla", required_argument, 0, 't'},
		{"finite", no_argument, &fitFinite, 1},
		{"contract", no_argument, 0, 'c'},
		{"resids", no_argument, 0, 'r'},
		{"threads", required_argument, 0, 'T'},
		{"no-prompt", no_argument, 0, 'n'},
        {"model", required_argument, 0, 'M'},
		{0, 0, 0, 0}
	};
    const char* short_options = "hvm:o:f:b:s:p:S:a:t:FT:crnM:i:j:";

	// Deal with these options in first pass to ensure the correct model is selected
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
            case 'M': {
                string choose_model(optarg);
                if (choose_model == "1") { model = make_shared<SCD>(); }
                else if (choose_model == "2") { model = make_shared<MCD2>(); }
                else if (choose_model == "2nex") { model = make_shared<MCD2_NoEx>(); }
                else if (choose_model == "3") { model = make_shared<MCD3>(); }
                else if (choose_model == "3nex") { model = make_shared<MCD3_NoEx>(); }
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
				mask = QI::ReadImageF::New();
				mask->SetFileName(optarg);
				break;
			case 'o':
				outPrefix = optarg;
				if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'f':
				if (verbose) cout << "Reading f0 file: " << optarg << endl;
				f0 = QI::ReadImageF::New();
				f0->SetFileName(optarg);
				break;
			case 'b':
				if (verbose) cout << "Reading B1 file: " << optarg << endl;
				B1 = QI::ReadImageF::New();
				B1->SetFileName(optarg);
				break;
			case 's': start_slice = atoi(optarg); break;
			case 'p': stop_slice = atoi(optarg); break;
			case 'S': {
				string mode(optarg);
				if (mode == "MEAN") {
					if (verbose) cout << "Mean scaling selected." << endl;
					model->setScaleToMean(true);
                    PDScale = 1.;
					applySlices->SetScaleToMean(true);
				} else if (atoi(optarg) == 0) {
					if (verbose) cout << "Fit PD/M0 selected." << endl;
                    model->setScaleToMean(false);
                    PDScale = atof(optarg);
					applySlices->SetScaleToMean(false);
				} else {
					if (verbose) cout << "Fix PD value to " << atof(optarg) << endl;
					model->setScaleToMean(false);
                    PDScale = atof(optarg);
					applySlices->SetScaleToMean(false);
				}
			} break;
            case 'a':
                switch (*optarg) {
                case 'S': algo = Algos::SRC; break;
                case 'G': algo = Algos::GRC; break;
                case 'L': algo = Algos::LBFGSB; nargs = 2; break;
                default:
                    cerr << "Unknown algorithm type " << *optarg << endl;
                    return EXIT_FAILURE;
                    break;
                } break;
			case 'F': flipData = true; break;
			case 'T': itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg)); break; break;
			case 't':
				switch (*optarg) {
					case '3': tesla = FieldStrength::Three; break;
					case '7': tesla = FieldStrength::Seven; break;
					case 'u': tesla = FieldStrength::User; break;
					default:
						cerr << "Unknown boundaries type " << *optarg << endl;
						return EXIT_FAILURE;
						break;
                } break;
			case 'c': {
				if (prompt) cout << "Enter max contractions/samples per contraction/retained samples/expand fraction: " << flush;
				ArrayXi in = ArrayXi::Zero(3);
				QI::ReadArray(cin, in);
                //mcd->setRCPars(in[0], in[1], in[2]);
			} break;
			case 'r': all_residuals = true; break;
			case 'h':
				cout << usage << endl;
				return EXIT_SUCCESS;
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
			case 0: break; // Just a flag
			default:
				cout << "Unhandled option " << string(1, c) << endl;
				return EXIT_FAILURE;
		}
	}
    if ((argc - optind) != nargs) {
		cerr << usage << endl << "Incorrect number of arguments." << endl;
		return EXIT_FAILURE;
	} else if (prompt) {
		cout << "Starting qimcdespot" << endl;
		cout << "Run with -h switch to see usage" << endl;
	}
    Array2d f0Bandwidth;

	shared_ptr<SequenceGroup> sequences = make_shared<SequenceGroup>();
	// Build a Functor here so we can query number of parameters etc.
	if (verbose) cout << "Using " << model->Name() << " model." << endl;
	vector<QI::ReadTimeseriesF::Pointer> inFiles;
	vector<QI::TimeseriesToVectorF::Pointer> inData;
	vector<QI::ReorderF::Pointer> inOrder;
	parseInput(sequences, inFiles, inData, inOrder, f0Bandwidth, fitFinite, flipData, verbose, prompt);

    ArrayXXd bounds = model->Bounds(tesla, 0);
    ArrayXd start = model->Start(tesla, 1., 1.);
    if (tesla == FieldStrength::User) {
        ArrayXd temp;
        if (prompt) cout << "Enter lower bounds" << endl;
        QI::ReadArray(cin, temp);
        bounds.col(0) = temp;
        if (prompt) cout << "Enter upper bounds" << endl;
        QI::ReadArray(cin, temp);
        bounds.col(1) = temp;
        if (algo == Algos::LBFGSB) {
            if (prompt) cout << "Enter start point" << endl;
            QI::ReadArray(cin, start);
        }
    }
    bounds.row(model->nParameters() - 2) = f0Bandwidth;
    shared_ptr<MCDAlgo> mcd;
    switch (algo) {
    case Algos::SRC:
        mcd = make_shared<SRCAlgo>();
        mcd->setModel(model);
        mcd->setBounds(bounds);
        mcd->setSequence(sequences);
        applySlices->SetAlgorithm(mcd);
        break;
    case Algos::GRC: {
        shared_ptr<SRCAlgo> temp = make_shared<SRCAlgo>();
        temp->setGauss(true);
        mcd = temp;
        mcd->setModel(model);
        mcd->setBounds(bounds);
        mcd->setSequence(sequences);
        applySlices->SetAlgorithm(mcd);
    } break;
    case Algos::LBFGSB: {
        itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);
        shared_ptr<LBFGSBAlgo> temp = make_shared<LBFGSBAlgo>();
        temp->setStart(start);
        QI::ReadImageF::Pointer T1 = QI::ReadImageF::New();
        QI::ReadImageF::Pointer T2 = QI::ReadImageF::New();
        T1->SetFileName(argv[optind++]);
        T2->SetFileName(argv[optind++]);
        T1->Update();
        T2->Update();
        mcd = temp;
        mcd->setModel(model);
        mcd->setBounds(bounds);
        mcd->setSequence(sequences);
        mcd->setTesla(tesla);
        applySlices->SetAlgorithm(mcd);
        applySlices->SetConst(2, T1->GetOutput());
        applySlices->SetConst(3, T2->GetOutput());
    } break;
    }
    mcd->setPDScale(PDScale);
	applySlices->SetSlices(start_slice, stop_slice);
	for (int i = 0; i < inOrder.size(); i++) {
		applySlices->SetInput(i, inOrder.at(i)->GetOutput());
	}
	if (f0) {
		f0->Update();
		applySlices->SetConst(0, f0->GetOutput());
	}
	if (B1) {
		B1->Update();
		applySlices->SetConst(1, B1->GetOutput());
	}
	if (mask) {
		mask->Update();
		applySlices->SetMask(mask->GetOutput());
	}

    // Need this here so the bounds.txt file will have the correct prefix
    outPrefix = outPrefix + model->Name() + "_";
    if (verbose) {
        cout << *sequences;
        cout << "Bounds:" << endl <<  bounds.transpose() << endl;
        if (algo == Algos::LBFGSB)
            cout << "Start: " << endl << start.transpose() << endl;
        ofstream boundsFile(outPrefix + "bounds.txt");
        boundsFile << "Names: ";
        for (size_t p = 0; p < model->nParameters(); p++) {
            boundsFile << model->Names()[p] << "\t";
        }
        boundsFile << endl << "Bounds: " << endl << bounds.transpose() << endl;
        if (algo == Algos::LBFGSB)
            boundsFile << "Start: " << endl << start.transpose() << endl;
        boundsFile.close();
    }

	time_t startTime;
	if (verbose) {
		startTime = QI::printStartTime();
		auto monitor = QI::SliceMonitor<TMCDFilter>::New();
		applySlices->AddObserver(itk::IterationEvent(), monitor);
	}
	applySlices->Update();
	if (verbose) {
		QI::printElapsedTime(startTime);
		cout << "Writing results files." << endl;
	}
	for (int i = 0; i < model->nParameters(); i++) {
		QI::writeResult(applySlices->GetOutput(i), outPrefix + model->Names()[i] + QI::OutExt());
	}
	QI::writeResiduals(applySlices->GetResidOutput(), outPrefix, all_residuals);
    QI::writeResult<itk::Image<int, 3>>(applySlices->GetIterationsOutput(), outPrefix + "iterations" + QI::OutExt());
	return EXIT_SUCCESS;
}

