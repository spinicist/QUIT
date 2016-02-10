/*
 *  despot2_main.cpp
 *
 *  Created by Tobias Wood on 23/01/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <getopt.h>
#include <iostream>
#include <atomic>

#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>
 
#include "Filters/ImageToVectorFilter.h"
#include "Filters/ApplyAlgorithmFilter.h"
#include "Model.h"
#include "Sequence.h"
#include "Util.h"
#include "itkTimeProbe.h"

using namespace std;
using namespace Eigen;
using namespace QI;

//******************************************************************************
// Algorithm Subclasses
//******************************************************************************
class D2Algo : public Algorithm<double> {
protected:
	const shared_ptr<SCD> m_model = make_shared<SCD>();
	shared_ptr<SteadyState> m_sequence;
	size_t m_iterations = 15;
	bool m_elliptical = false;
	double m_thresh = -numeric_limits<double>::infinity();
	double m_lo = -numeric_limits<double>::infinity();
	double m_hi = numeric_limits<double>::infinity();

public:
	void setIterations(size_t n) { m_iterations = n; }
	void setSequence(shared_ptr<SteadyState> &s) { m_sequence = s; }
	void setElliptical(bool e) { m_elliptical = e; }
	void setThreshold(double t) { m_thresh = t; }
	void setClamp(double lo, double hi) { m_lo = lo; m_hi = hi; }
	size_t numInputs() const override { return m_sequence->count(); }
	size_t numConsts() const override { return 2; }  // T1, B1
	size_t numOutputs() const override { return 2; } // PD, T2
	size_t dataSize() const override { return m_sequence->size(); }

    virtual TArray defaultConsts() override {
		// T1, B1
		TArray def = TArray::Ones(2);
		return def;
	}
};

class D2LLS : public D2Algo {
public:
	virtual void apply(const TInput &data, const TArray &constants,
                       TArray &outputs, TArray &resids, TIterations &its) const override
	{
		const double TR = m_sequence->TR();
		const double T1 = constants[0];
		const double B1 = constants[1];
		const double E1 = exp(-TR / T1);
		double PD, T2, E2;
		const ArrayXd angles = (m_sequence->flip() * B1);

		VectorXd Y = data / angles.sin();
		MatrixXd X(Y.rows(), 2);
		X.col(0) = data / angles.tan();
		X.col(1).setOnes();
		VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
		if (m_elliptical) {
			T2 = 2. * TR / log((b[0]*E1 - 1.) / (b[0] - E1));
			E2 = exp(-TR / T2);
			PD = b[1] * (1. - E1*E2*E2) / (sqrt(E2) * (1. - E1));
		} else {
			T2 = TR / log((b[0]*E1 - 1.)/(b[0] - E1));
			E2 = exp(-TR / T2);
			PD = b[1] * (1. - E1*E2) / (1. - E1);
		}
		VectorXd p(5); p << PD, T1, T2, 0, B1;
		ArrayXd theory = m_sequence->signal(m_model, p).abs();
		resids = data.array() - theory;
		if (PD < m_thresh)
			PD = T2 = 0.;
		T2 = clamp(T2, m_lo, m_hi);
		outputs[0] = PD;
		outputs[1] = T2;
        its = 1;
	}
};

class D2WLLS : public D2Algo {
public:
	virtual void apply(const TInput &data, const TArray &constants,
                       TArray &outputs, TArray &resids, TIterations &its) const override
	{
		const double TR = m_sequence->TR();
		const double T1 = constants[0];
		const double B1 = constants[1];
		const double E1 = exp(-TR / T1);
		double PD, T2, E2;
		const ArrayXd angles = (m_sequence->flip() * B1);

		VectorXd Y = data / angles.sin();
		MatrixXd X(Y.rows(), 2);
		X.col(0) = data / angles.tan();
		X.col(1).setOnes();
		VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
		if (m_elliptical) {
			T2 = 2. * TR / log((b[0]*E1 - 1.) / (b[0] - E1));
			E2 = exp(-TR / T2);
			PD = b[1] * (1. - E1*E2*E2) / (sqrt(E2) * (1. - E1));
		} else {
			T2 = TR / log((b[0]*E1 - 1.)/(b[0] - E1));
			E2 = exp(-TR / T2);
			PD = b[1] * (1. - E1*E2) / (1. - E1);
		}
		VectorXd W(m_sequence->size());
		for (size_t n = 0; n < m_iterations; n++) {
			if (m_elliptical) {
				W = ((1. - E1*E2) * angles.sin() / (1. - E1*E2*E2 - (E1 - E2*E2)*angles.cos())).square();
			} else {
				W = ((1. - E1*E2) * angles.sin() / (1. - E1*E2 - (E1 - E2)*angles.cos())).square();
			}
			b = (X.transpose() * W.asDiagonal() * X).partialPivLu().solve(X.transpose() * W.asDiagonal() * Y);
			if (m_elliptical) {
				T2 = 2. * TR / log((b[0]*E1 - 1.) / (b[0] - E1));
				E2 = exp(-TR / T2);
				PD = b[1] * (1. - E1*E2*E2) / (sqrt(E2) * (1. - E1));
			} else {
				T2 = TR / log((b[0]*E1 - 1.)/(b[0] - E1));
				E2 = exp(-TR / T2);
				PD = b[1] * (1. - E1*E2) / (1. - E1);
			}
		}
		VectorXd p(5); p << PD, T1, T2, 0, B1;
		ArrayXd theory = m_sequence->signal(m_model, p).abs();
		resids = data.array() - theory;
		if (PD < m_thresh)
			PD = T2 = 0.;
		T2 = clamp(T2, m_lo, m_hi);
		outputs[0] = PD;
		outputs[1] = T2;
        its = m_iterations;
	}
};

//******************************************************************************
// T2 Only Functor
//******************************************************************************
class D2Functor : public DenseFunctor<double> {
	public:
		const shared_ptr<SequenceBase> m_sequence;
		const double m_T1, m_B1;
		const shared_ptr<SCD> m_model = make_shared<SCD>();
		const ArrayXd m_data;

		D2Functor(const double T1, const shared_ptr<SequenceBase> s, const ArrayXd &d, const double B1, const bool fitComplex, const bool debug = false) :
			DenseFunctor<double>(3, s->size()),
			m_sequence(s), m_data(d),
			m_T1(T1), m_B1(B1)
		{
			assert(static_cast<size_t>(m_data.rows()) == values());
		}

		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
			eigen_assert(diffs.size() == values());
			ArrayXd fullparams(5);
			fullparams << params(0), m_T1, params(1), params(2), m_B1;
			ArrayXcd s = m_sequence->signal(m_model, fullparams);
			diffs = s.abs() - m_data;
			return 0;
		}
};

class D2NLLS : public D2Algo {
public:
	virtual void apply(const TInput &data, const TArray &inputs,
                       TArray &outputs, TArray &resids, TIterations &its) const override
	{
		double T1 = inputs[0];
		double B1 = inputs[1];
		D2Functor f(T1, m_sequence, data, B1, false, false);
		NumericalDiff<D2Functor> nDiff(f);
		LevenbergMarquardt<NumericalDiff<D2Functor>> lm(nDiff);
		lm.setMaxfev(m_iterations * (m_sequence->size() + 1));
		VectorXd p(2); p << data.array().maxCoeff() * 5., 0.1;
		lm.minimize(p);
		outputs = p;
		if (outputs[0] < m_thresh)
			outputs.setZero();
		outputs[1] = clamp(outputs[1], m_lo, m_hi);
		VectorXd fullp(5); fullp << outputs[0], T1, outputs[1], 0, B1; // Assume on-resonance
		ArrayXd theory = m_sequence->signal(m_model, fullp).abs(); // Sequence will already be elliptical if necessary
		resids = data.array() - theory;
		if (outputs[0] < m_thresh)
			outputs.setZero();
		outputs[1] = clamp(outputs[1], m_lo, m_hi);
        its = lm.iterations();
	}
};

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: qidespot2 [options] T1_map ssfp_file\n\
\n\
Options:\n\
	--help, -h        : Print this message.\n\
	--verbose, -v     : Print slice processing times.\n\
	--no-prompt, -n   : Suppress input prompts.\n\
	--mask, -m file   : Mask input with specified file.\n\
	--out, -o path    : Add a prefix to the output filenames.\n\
	--B1 file         : B1 Map file.\n\
	--thresh, -t n    : Threshold maps at PD < n\n\
	--clamp, -c n     : Clamp T2 between 0 and n\n\
	--algo, -a l      : LLS algorithm (default)\n\
	           w      : WLLS algorithm\n\
	           n      : NLLS (Levenberg-Marquardt)\n\
	--its, -i N       : Max iterations for WLLS / NLLS (default 10)\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n\
	--elliptical, -e  : Input is band-free elliptical data.\n"
};

int verbose = false, prompt = true, elliptical = false, all_residuals = false, num_threads = 4;
string outPrefix;
const struct option long_opts[] = {
	{"B1", required_argument, 0, 'b'},
	{"elliptical", no_argument, 0, 'e'},
	{"help", no_argument, 0, 'h'},
	{"mask", required_argument, 0, 'm'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"thresh", required_argument, 0, 't'},
	{"clamp", required_argument, 0, 'c'},
	{"algo", required_argument, 0, 'a'},
	{"its", required_argument, 0, 'i'},
	{"threads", required_argument, 0, 'T'},
	{"resids", no_argument, 0, 'r'},
	{0, 0, 0, 0}
};
const char *short_opts = "hm:o:b:t:c:vna:i:T:er";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();
	ImageReaderF::Pointer mask = ITK_NULLPTR;
	ImageReaderF::Pointer B1   = ITK_NULLPTR;
	shared_ptr<D2Algo> algo = make_shared<D2LLS>();

	// Do first pass to get the algorithm type, then do everything else
	int indexptr = 0;
	char c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'a':
			switch (*optarg) {
				case 'l': algo = make_shared<D2LLS>();   if (verbose) cout << "LLS algorithm selected." << endl; break;
				case 'w': algo = make_shared<D2WLLS>();  if (verbose) cout << "WLLS algorithm selected." << endl; break;
				case 'n': algo = make_shared<D2NLLS>(); if (verbose) cout << "NLLS algorithm selected." << endl; break;
				default: throw(runtime_error(string("Unknown algorithm type ") + optarg)); break;
			} break;
			default: break;
		}
	}

	optind = 1;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
			case 'v': case 'n': case 'a': break; // Already handled
			case 'o':
				outPrefix = optarg;
				if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'm':
				if (verbose) cout << "Reading mask file " << optarg << endl;
				mask = ImageReaderF::New();
				mask->SetFileName(optarg);
				break;
			case 'b':
				if (verbose) cout << "Reading B1 file: " << optarg << endl;
				B1 = ImageReaderF::New();
				B1->SetFileName(optarg);
				break;
			case 't': algo->setThreshold(atof(optarg)); break;
			case 'c': algo->setClamp(0, atof(optarg)); break;
			case 'i': algo->setIterations(atoi(optarg));	break;
			case 'e': elliptical = true; break;
			case 'r': all_residuals = true; break;
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
	//if (verbose) cout << version << endl << credit_shared << endl;
	if ((argc - optind) != 2) {
		cout << "Wrong number of arguments. Need a T1 map and SSFP file." << endl;
		cout << usage << endl;
		return EXIT_FAILURE;
	}
	if (verbose) cout << "Reading T1 Map from: " << argv[optind] << endl;
	auto T1File = ImageReaderF::New();
	T1File->SetFileName(argv[optind++]);

	if (verbose) cout << "Opening SSFP file: " << argv[optind] << endl;
	auto ssfp4D = TimeseriesReaderF::New();
	ssfp4D->SetFileName(argv[optind++]);
	auto ssfp3D = itk::ImageToVectorFilter<TimeseriesF>::New();
	ssfp3D->SetInput(ssfp4D->GetOutput());

	shared_ptr<SteadyState> ssfp;
	if (elliptical) {
		ssfp = make_shared<SSFPEllipse>(prompt);
	} else {
		ssfp = make_shared<SSFPSimple>(prompt);
	}
	if (verbose) cout << *ssfp << endl;

    auto DESPOT2 = itk::ApplyAlgorithmFilter<D2Algo>::New();
	algo->setSequence(ssfp);
	algo->setElliptical(elliptical);
	DESPOT2->SetAlgorithm(algo);
    DESPOT2->SetPoolsize(num_threads);
	DESPOT2->SetInput(0, ssfp3D->GetOutput());
	DESPOT2->SetConst(0, T1File->GetOutput());
	if (B1)
		DESPOT2->SetConst(1, B1->GetOutput());
	if (mask)
		DESPOT2->SetMask(mask->GetOutput());

    itk::TimeProbe clock;
	if (verbose) {
		cout << "DESPOT2 setup complete. Processing." << endl;
		auto monitor = QI::GenericMonitor::New();
		DESPOT2->AddObserver(itk::ProgressEvent(), monitor);
        clock.Start();
	}
	DESPOT2->Update();
	if (verbose) {
        clock.Stop();
        clock.Stop();
        cout << "Elapsed time was " << clock.GetTotal() << "s" << endl;
		cout << "Writing results." << endl;
	}
	outPrefix = outPrefix + "D2_";
    WriteImage(DESPOT2->GetOutput(0), outPrefix + "PD.nii");
	WriteImage(DESPOT2->GetOutput(1), outPrefix + "T2.nii");
	WriteResiduals(DESPOT2->GetResidOutput(), outPrefix, all_residuals, DESPOT2->GetOutput(0));

	if (verbose) cout << "All done." << endl;
	return EXIT_SUCCESS;
}
