/*
 *  multiecho.cpp
 *
 *  Created by Tobias Wood on 27/01/2015.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <time.h>
#include <getopt.h>
#include <iostream>
#include <atomic>
#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "Types.h"
#include "Util.h"
#include "Sequence.h"
#include "Filters/ApplyAlgorithmFilter.h"

using namespace std;
using namespace Eigen;

const string usage {
"Usage is: multiecho [options] input_file \n\
\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--no-prompt, -n   : Suppress input prompts\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--mask, -m file   : Mask input with specified file\n\
	--star, -S        : Data is T2*, not T2\n\
	--thresh, -t n    : Threshold maps at PD < n\n\
	--clamp, -c n     : Clamp T2 between 0 and n\n\
	--algo, -a l      : LLS algorithm (default)\n\
	           a      : ARLO algorithm\n\
	           n      : Non-linear (Levenberg-Marquardt)\n\
	--its, -i N       : Max iterations for non-linear (default 10)\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

static int NE = 0, nIterations = 10;
static bool verbose = false, prompt = true, all_residuals = false, weightedSum = false;
static string outPrefix, suffix;
static double thresh = -numeric_limits<double>::infinity();
static double clamp_lo = -numeric_limits<double>::infinity(), clamp_hi = numeric_limits<double>::infinity();
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"star", no_argument, 0, 'S'},
	{"thresh", required_argument, 0, 't'},
	{"clamp", required_argument, 0, 'c'},
	{"algo", required_argument, 0, 'a'},
	{"threads", required_argument, 0, 'T'},
	{"resids", no_argument, 0, 'r'},
	{0, 0, 0, 0}
};
static const char *short_opts = "hvnm:Se:o:b:t:c:a:T:r";

/*
 * Base class for the 3 different algorithms
 */
class RelaxAlgo : public Algorithm<double> {
private:
	const shared_ptr<SCD> m_model = make_shared<SCD>();
protected:
	shared_ptr<MultiEcho> m_sequence;
    double m_clampLo = -numeric_limits<double>::infinity();
    double m_clampHi = numeric_limits<double>::infinity();
    double m_thresh = -numeric_limits<double>::infinity();

    void clamp_and_treshold(const TInput &data, TArray &outputs, TArray &resids,
                            const double PD, const double T2) const {
        if (PD > m_thresh) {
            outputs[0] = PD;
            outputs[1] = clamp(T2, m_clampLo, m_clampHi);
            ArrayXcd theory = One_MultiEcho(m_sequence->m_TE, PD, T2);
            resids = data.array() - theory.abs();
        } else {
            outputs.setZero();
            resids.setZero();
        }
    }

public:
	void setSequence(shared_ptr<MultiEcho> &s) { m_sequence = s; }
    void setClamp(double lo, double hi) { m_clampLo = lo; m_clampHi = hi; }
    void setThresh(double t) { m_thresh = t; }
    size_t numInputs() const override { return m_sequence->count(); }
	size_t numConsts() const override { return 1; }
	size_t numOutputs() const override { return 2; }
	size_t dataSize() const override { return m_sequence->size(); }

	virtual TArray defaultConsts() {
		// B1
		TArray def = TArray::Ones(2);
		return def;
	}
};

class LogLinAlgo: public RelaxAlgo {
public:
	virtual void apply(const TInput &data, const TArray &inputs,
                       TArray &outputs, TArray &resids, TIterations &its) const override
	{
			// Set up echo times array
		MatrixXd X(m_sequence->size(), 2);
		X.col(0) = m_sequence->m_TE;
		X.col(1).setOnes();
		VectorXd Y = data.array().log();
		VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        double PD = exp(b[1]);
        double T2 = -1 / b[0];
        clamp_and_treshold(data, outputs, resids, PD, T2);
        its = 1;
	}
};

class ARLOAlgo : public RelaxAlgo {
public:
    virtual void apply(const TInput &data, const TArray &inputs,
                       TArray &outputs, TArray &resids, TIterations &its) const override
    {
        const double dTE_3 = (m_sequence->m_ESP / 3);
        double si2sum = 0, sidisum = 0;
        for (int i = 0; i < m_sequence->size() - 2; i++) {
            const double si = dTE_3 * (data(i) + 4*data(i+1) + data(i+2));
            const double di = data(i) - data(i+2);
            si2sum += si*si;
            sidisum += si*di;
        }
        double T2 = (si2sum + dTE_3*sidisum) / (dTE_3*si2sum + sidisum);
        double PD = (data.array() / exp(-m_sequence->m_TE / T2)).mean();
        clamp_and_treshold(data, outputs, resids, PD, T2);
        its = 1;
    }
};

class RelaxFunctor : public DenseFunctor<double> {
	protected:
		const shared_ptr<SequenceBase> m_sequence;
		const ArrayXd m_data;
		const shared_ptr<SCD> m_model = make_shared<SCD>();

	public:
		RelaxFunctor(shared_ptr<SequenceBase> cs, const ArrayXd &data) :
			DenseFunctor<double>(2, cs->size()),
			m_sequence(cs), m_data(data)
		{
			assert(static_cast<size_t>(m_data.rows()) == values());
		}

		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
			eigen_assert(diffs.size() == values());
			VectorXd fullp(5);
			fullp << params(0), 0, params(1), 0, 1.0; // Fix B1 to 1.0 for now
			ArrayXcd s = m_sequence->signal(m_model, fullp);
			diffs = s.abs() - m_data;
			return 0;
		}
};

class NonLinAlgo : public RelaxAlgo {
private:
	size_t m_iterations = 5;
public:
	void setIterations(size_t n) { m_iterations = n; }

	virtual void apply(const TInput &data, const TArray &inputs,
                       TArray &outputs, TArray &resids, TIterations &its) const override
	{
		RelaxFunctor f(m_sequence, data);
		NumericalDiff<RelaxFunctor> nDiff(f);
		LevenbergMarquardt<NumericalDiff<RelaxFunctor>> lm(nDiff);
		lm.setMaxfev(nIterations * (m_sequence->size() + 1));
		// Just PD & T2 for now
		// Basic guess of T2=50ms
		VectorXd p(2); p << data(0), 0.05;
		lm.minimize(p);
        clamp_and_treshold(data, outputs, resids, p[0], p[1]);
        its = lm.iterations();
	}
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();
	QI::ReadImageF::Pointer mask, B1, f0 = ITK_NULLPTR;
	shared_ptr<RelaxAlgo> algo = make_shared<LogLinAlgo>();
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				mask = QI::ReadImageF::New();
				mask->SetFileName(optarg);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'S': suffix = "star"; break;
			case 't': thresh = atof(optarg); break;
			case 'c':
				clamp_lo = 0;
				clamp_hi = atof(optarg);
				break;
			case 'a':
				switch (*optarg) {
					case 'l': algo = make_shared<LogLinAlgo>(); if (verbose) cout << "LogLin algorithm selected." << endl; break;
					case 'a': algo = make_shared<ARLOAlgo>(); if (verbose) cout << "ARLO algorithm selected." << endl; break;
					case 'n': algo = make_shared<NonLinAlgo>(); if (verbose) cout << "Non-linear algorithm (Levenberg Marquardt) selected." << endl; break;
					default:
						cout << "Unknown algorithm type " << optarg << endl;
						return EXIT_FAILURE;
						break;
				} break;
			case 'T': itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg)); break;
			case 'r': all_residuals = true; break;
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
	if ((argc - optind) != 1) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		return EXIT_FAILURE;
	}
	// Gather input data
	cout << "Opening input file: " << argv[optind] << endl;
	auto inputFile = QI::ReadTimeseriesF::New();
	inputFile->SetFileName(argv[optind]);
	auto inputData = QI::TimeseriesToVectorF::New();
	inputData->SetInput(inputFile->GetOutput());
	auto multiecho = make_shared<MultiEcho>(prompt);
	if (verbose) {
		cout << "Ouput prefix will be: " << outPrefix << endl;
		cout << "Clamp: " << clamp_lo << " " << clamp_hi << endl;
		cout << "Thresh: " << thresh << endl;
	}

	algo->setSequence(multiecho);
    auto apply = itk::ApplyAlgorithmFilter<RelaxAlgo>::New();
	apply->SetAlgorithm(algo);
	apply->SetInput(0, inputData->GetOutput());
	if (mask)
		apply->SetMask(mask->GetOutput());

	time_t startTime;
	if (verbose) {
		startTime = QI::printStartTime();
		auto monitor = QI::GenericMonitor::New();
		apply->AddObserver(itk::ProgressEvent(), monitor);
	}
	apply->Update();
	if (verbose) {
		QI::printElapsedTime(startTime);
		cout << "Writing results." << endl;
	}
	outPrefix = outPrefix + "ME_";
    QI::writeResult(apply->GetOutput(0), outPrefix + "PD" + suffix + QI::OutExt());
    QI::writeResult(apply->GetOutput(1), outPrefix + "T2" + suffix + QI::OutExt());
    QI::writeResiduals(apply->GetResidOutput(), outPrefix, all_residuals);

	return EXIT_SUCCESS;
}

