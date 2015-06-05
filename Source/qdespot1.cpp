/*
 *  qdespot1.cpp - Part of QUantitative Imaging Tools
 *
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

#include "itkImageFileReader.h"

#include "Filters/ImageToVectorFilter.h"
#include "Filters/ApplyAlgorithmFilter.h"
#include "Model.h"
#include "Sequence.h"
#include "Util.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
// Algorithm Subclasses
//******************************************************************************
class D1Algo : public Algorithm<double> {
protected:
	const shared_ptr<Model> m_model = make_shared<SCD>();
	shared_ptr<SPGRSimple> m_sequence;
	size_t m_iterations = 15;

public:
	void setIterations(size_t n) { m_iterations = n; }
	void setSequence(shared_ptr<SPGRSimple> &s) { m_sequence = s; }
	size_t numInputs() const override { return m_sequence->count(); }
	size_t numConsts() const override { return 1; }
	size_t numOutputs() const override { return 2; }
	size_t dataSize() const override { return m_sequence->size(); }

	virtual VectorXd defaultConsts() {
		// B1
		VectorXd def = VectorXd::Ones(1);
		return def;
	}
};

class D1LLS : public D1Algo {
public:
	virtual void apply(const VectorXd &data, const VectorXd &inputs,
	                   VectorXd &outputs, ArrayXd &resids) const override
	{
		double B1 = inputs[0];
		ArrayXd flip = m_sequence->flip() * B1;
		VectorXd Y = data.array() / flip.sin();
		MatrixXd X(Y.rows(), 2);
		X.col(0) = data.array() / flip.tan();
		X.col(1).setOnes();
		VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
		outputs[1] = -m_sequence->TR() / log(b[0]);
		outputs[0] = b[1] / (1. - b[0]);
		ArrayXd theory = One_SPGR(m_sequence->flip(), m_sequence->TR(), outputs[0], outputs[1], B1).array().abs();
		resids = (data.array() - theory);
	}
};

class D1WLLS : public D1Algo {
public:
	virtual void apply(const VectorXd &data, const VectorXd &inputs,
	                   VectorXd &outputs, ArrayXd &resids) const override
	{
		double B1 = inputs[0];
		ArrayXd flip = m_sequence->flip() * B1;
		VectorXd Y = data.array() / flip.sin();
		MatrixXd X(Y.rows(), 2);
		X.col(0) = data.array() / flip.tan();
		X.col(1).setOnes();
		VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
		outputs[1] = -m_sequence->TR() / log(b[0]);
		outputs[0] = b[1] / (1. - b[0]);
		for (size_t n = 0; n < m_iterations; n++) {
			MatrixXd W = (flip.sin() / (1. - (exp(-m_sequence->TR()/outputs[1])*flip.cos()))).square();
			b = (X.transpose() * W.asDiagonal() * X).partialPivLu().solve(X.transpose() * W.asDiagonal() * Y);
			outputs[1] = -m_sequence->TR() / log(b[0]);
			outputs[0] = b[1] / (1. - b[0]);
		}
		ArrayXd theory = One_SPGR(m_sequence->flip(), m_sequence->TR(), outputs[0], outputs[1], B1).array().abs();
		resids = (data.array() - theory);
	}
};

// T1 only Functor
class T1Functor : public DenseFunctor<double> {
	protected:
		const shared_ptr<SequenceBase> m_sequence;
		const ArrayXd m_data;
		const double m_B1;
		const shared_ptr<SCD> m_model;

	public:
		T1Functor(const shared_ptr<SequenceBase> cs, const ArrayXd &data, const double B1) :
			DenseFunctor<double>(2, cs->size()),
			m_sequence(cs), m_data(data), m_B1(B1)
		{
			assert(static_cast<size_t>(m_data.rows()) == values());
		}

		int operator()(const Ref<VectorXd> &p, Ref<ArrayXd> diffs) const {
			eigen_assert(diffs.size() == values());
			ArrayXd s = One_SPGR(m_sequence->flip(), m_sequence->TR(), p[0], p[1], m_B1).array().abs();
			diffs = s.abs() - m_data;
			return 0;
		}
};

class D1NLLS : public D1Algo {
public:
	virtual void apply(const VectorXd &data, const VectorXd &inputs,
	                   VectorXd &outputs, ArrayXd &resids) const override
	{
		double B1 = inputs[0];
		T1Functor f(m_sequence, data, B1);
		NumericalDiff<T1Functor> nDiff(f);
		LevenbergMarquardt<NumericalDiff<T1Functor>> lm(nDiff);
		// PD & T1 - Initial guess of 1s
		outputs << data.array().abs().maxCoeff() * 10., 1.;
		lm.setMaxfev(m_iterations * (m_sequence->size() + 1));
		lm.minimize(outputs);
		ArrayXd theory = One_SPGR(m_sequence->flip(), m_sequence->TR(), outputs[0], outputs[1], B1).array().abs();
		resids = data.array() - theory;
	}
};

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: despot1 [options] spgr_input \n\
\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--no-prompt, -n   : Suppress input prompts\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--mask, -m file   : Mask input with specified file\n\
	--B1, -b file     : B1 Map file (ratio)\n\
	--algo, -a l      : LLS algorithm (default)\n\
	           w      : WLLS algorithm\n\
	           n      : NLLS (Levenberg-Marquardt)\n\
	--its, -i N       : Max iterations for WLLS/NLLS (default 15)\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

static bool verbose = false, prompt = true, all_residuals = false;
static size_t nIterations = 4;
static string outPrefix;
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"B1", required_argument, 0, 'b'},
	{"algo", required_argument, 0, 'a'},
	{"its", required_argument, 0, 'i'},
	{"resids", no_argument, 0, 'r'},
	{"threads", required_argument, 0, 'T'},
	{0, 0, 0, 0}
};
static const char *short_opts = "hvnm:o:b:a:i:rT:";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	try { // To fix uncaught exceptions on Mac
	//cout << version << endl << credit_shared << endl;
	Eigen::initParallel();

	QI::ReadImageF::Pointer mask = ITK_NULLPTR;
	QI::ReadImageF::Pointer B1   = ITK_NULLPTR;

	shared_ptr<D1Algo> algo = make_shared<D1LLS>();
	size_t its = 15;
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'a':
				switch (*optarg) {
					case 'l': algo = make_shared<D1LLS>();  if (verbose) cout << "LLS algorithm selected." << endl; break;
					case 'w': algo = make_shared<D1WLLS>(); if (verbose) cout << "WLLS algorithm selected." << endl; break;
					case 'n': algo = make_shared<D1NLLS>(); if (verbose) cout << "NLLS algorithm selected." << endl; break;
					default:
						cerr << "Unknown algorithm type " << optarg << endl;
						return EXIT_FAILURE;
						break;
				} break;
			case 'm':
				if (verbose) cout << "Opening mask file " << optarg << endl;
				mask = QI::ReadImageF::New();
				mask->SetFileName(optarg);
				break;
			case 'o':
				outPrefix = optarg;
				if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'b':
				if (verbose) cout << "Opening B1 file: " << optarg << endl;
				B1 = QI::ReadImageF::New();
				B1->SetFileName(optarg);
				break;
			case 'i': its = (atoi(optarg)); break;
			case 'r': all_residuals = true; break;
			case 'T': itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg)); break;
			case 'h':
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 1) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		return EXIT_FAILURE;
	}

	string inputFilename = argv[optind++];
	if (verbose) cout << "Opening SPGR file: " << inputFilename << endl;
	auto input = QI::ReadTimeseriesF::New();
	auto convert = QI::TimeseriesToVectorF::New();
	input->SetFileName(inputFilename);
	convert->SetInput(input->GetOutput());

	shared_ptr<SPGRSimple> spgrSequence = make_shared<SPGRSimple>(prompt);
	if (verbose) cout << *spgrSequence;
	algo->setSequence(spgrSequence);
	algo->setIterations(its);
	auto apply = itk::ApplyAlgorithmFilter<QI::VectorImageF, D1Algo>::New();
	apply->SetAlgorithm(algo);
	apply->SetDataInput(0, convert->GetOutput());
	if (mask)
		apply->SetMask(mask->GetOutput());
	if (B1)
		apply->SetConstInput(0, B1->GetOutput());
	if (verbose) cout << "Processing" << endl;
	apply->Update();
	if (verbose) cout << "Writing results." << endl;
	outPrefix = outPrefix + "D1_";

	QI::writeResult(apply->GetOutput(0), outPrefix + "PD.nii");
	QI::writeResult(apply->GetOutput(1), outPrefix + "T1.nii");
	QI::writeResiduals(apply->GetResidOutput(), outPrefix, all_residuals);

	if (verbose) cout << "Finished." << endl;
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
