/*
 *  despot1hifi.cpp
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright (c) 2011-2013 Tobias Wood.
 *
 *  Based on code by Sean Deoni
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

#include "itkImageFileReader.h"

#include "Sequence.h"
#include "Util.h"
#include "Filters/ImageToVectorFilter.h"
#include "Filters/ApplyAlgorithmFilter.h"

using namespace std;
using namespace Eigen;
using namespace QI;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: despot1hifi [options] spgr_input ir-spgr_input\n\
\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--no-prompt, -n   : Suppress input prompts\n\
	--mprage, -M      : Use a generic MP-RAGE sequence, not GE IR-SPGR\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--mask, -m file   : Mask input with specified file\n\
	--thresh, -t n    : Threshold maps at PD < n\n\
	--clamp, -c n     : Clamp T1 between 0 and n\n\
	--its, -i N       : Max iterations for NLLS (default 4)\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

static bool verbose = false, prompt = true, IR = true, all_residuals = false;
static size_t nIterations = 4;
static string outPrefix;
static double thresh = -numeric_limits<double>::infinity();
static double clamp_lo = -numeric_limits<double>::infinity(), clamp_hi = numeric_limits<double>::infinity();
static const struct option long_opts[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"mprage", no_argument, 0, 'M'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"thresh", required_argument, 0, 't'},
	{"clamp", required_argument, 0, 'c'},
	{"its", required_argument, 0, 'i'},
	{"resids", no_argument, 0, 'r'},
	{"threads", required_argument, 0, 'T'},
	{0, 0, 0, 0}
};
static const char *short_opts = "hvnMm:o:t:c:s:p:i:rT:";

// HIFI Algorithm - includes optimising B1
class HIFIFunctor : public DenseFunctor<double> {
	protected:
		const shared_ptr<SequenceBase> m_sequence;
		const ArrayXd m_data;
		const shared_ptr<SCD> m_model = make_shared<SCD>();

	public:
		HIFIFunctor(const shared_ptr<SequenceBase> cs, const ArrayXd &data) :
			DenseFunctor<double>(3, cs->size()),
			m_sequence(cs), m_data(data)
		{
			assert(static_cast<size_t>(m_data.rows()) == values());
		}

		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
			eigen_assert(diffs.size() == values());
			VectorXd fullpar(5); // Add T2 and f0
			fullpar << params(0), params(1), 0, 0, params(2);
			ArrayXcd s = m_sequence->signal(m_model, fullpar);
			diffs = s.abs() - m_data;
			return 0;
		}
};

class HIFIAlgo : public Algorithm<double> {
	private:
		shared_ptr<SequenceGroup> m_sequence;
		size_t m_iterations = 15; // From tests this seems to be a sensible maximum number
		double m_thresh = -numeric_limits<double>::infinity();
		double m_lo = -numeric_limits<double>::infinity();
		double m_hi = numeric_limits<double>::infinity();
	public:
		void setSequence(shared_ptr<SequenceGroup> & s) { m_sequence = s; }
		void setIterations(size_t n) { m_iterations = n; }
		void setThreshold(double t) { m_thresh = t; }
		void setClamp(double lo, double hi) { m_lo = lo; m_hi = hi; }
		size_t numInputs() const override  { return m_sequence->count(); }
		size_t numConsts() const override  { return 0; }
		size_t numOutputs() const override { return 3; }
		size_t dataSize() const override   { return m_sequence->size(); }

		virtual TArray defaultConsts() {
			// No constants for HIFI
			TArray def = TArray::Zero(0);
			return def;
		}

		virtual void apply(const TInput &data,
		                   const TArray &, //No inputs, remove name to silence compiler warning
		                   TArray &outputs, TArray &resids) const override
		{
			HIFIFunctor f(m_sequence, data);
			NumericalDiff<HIFIFunctor> nDiff(f);
			LevenbergMarquardt<NumericalDiff<HIFIFunctor>> lm(nDiff);
			// LevenbergMarquardt does not currently have a good interface, have to do things in steps
			lm.setMaxfev(m_iterations * (m_sequence->size() + 1));
			VectorXd op(3); op << data.array().abs().maxCoeff() * 10., 1., 1.; // Initial guess
			lm.minimize(op);
			outputs = op;
			// PD, T1, B1
			VectorXd pfull(5); pfull << outputs[0], outputs[1], 0, 0, outputs[2]; // Build full parameter vector
			auto model = make_shared<SCD>();
			ArrayXd theory = m_sequence->signal(model, pfull).abs();
			resids = (data.array() - theory);
			if (outputs[0] < m_thresh)
				outputs.setZero();
			outputs[1] = clamp(outputs[1], m_lo, m_hi);
		}
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();
	ReadImageF::Pointer mask = ITK_NULLPTR;
	auto hifi = make_shared<HIFIAlgo>();
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'M': IR = false; break;
			case 'm':
				if (verbose) cout << "Opening mask file: " << optarg << endl;
				mask = ReadImageF::New();
				mask->SetFileName(optarg);
				break;
			case 'o':
				outPrefix = optarg;
				if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 't': hifi->setThreshold(atof(optarg)); break;
			case 'c': hifi->setClamp(0, atof(optarg)); break;
			case 'i': hifi->setIterations(atoi(optarg)); break;
			case 'r': all_residuals = true; break;
			case 'T':
				itk::MultiThreader::SetGlobalDefaultNumberOfThreads(atoi(optarg));
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
	if ((argc - optind) != 2) {
		cerr << "Incorrect number of arguments." << endl;
		cout << usage << endl;
		return EXIT_FAILURE;
	}
	
	if (verbose) cout << "Opening SPGR file: " << argv[optind] << endl;
	auto spgrFile = ReadTimeseriesF::New();
	auto spgrImg = TimeseriesToVectorF::New();
	spgrFile->SetFileName(argv[optind++]);
	spgrImg->SetInput(spgrFile->GetOutput());
	auto spgrSequence = make_shared<SPGRSimple>(prompt);
	if (verbose) cout << "Opening IR-SPGR file: " << argv[optind] << endl;
	auto irFile = ReadTimeseriesF::New();
	auto irImg = TimeseriesToVectorF::New();
	irFile->SetFileName(argv[optind++]);
	irImg->SetInput(irFile->GetOutput());
	shared_ptr<SteadyState> irSequence;
	if (IR) {
		irSequence = make_shared<IRSPGR>(prompt);
	} else {
		irSequence = make_shared<MPRAGE>(prompt);
	}
	auto combined = make_shared<SequenceGroup>();
	combined->addSequence(spgrSequence);
	combined->addSequence(irSequence);
	if (verbose) cout << *combined << endl;
	auto apply = itk::ApplyAlgorithmFilter<QI::VectorImageF, HIFIAlgo>::New();
	hifi->setSequence(combined);
	apply->SetAlgorithm(hifi);
	apply->SetInput(0, spgrImg->GetOutput());
	apply->SetInput(1, irImg->GetOutput());
	if (mask)
		apply->SetMask(mask->GetOutput());
	if (verbose) {
		cout << "Processing..." << endl;
		auto monitor = QI::GenericMonitor::New();
		apply->AddObserver(itk::ProgressEvent(), monitor);
	}
	apply->Update();
	if (verbose) cout << "Writing results." << endl;
	outPrefix = outPrefix + "HIFI_";

	QI::writeResult(apply->GetOutput(0), outPrefix + "PD.nii");
	QI::writeResult(apply->GetOutput(1), outPrefix + "T1.nii");
	QI::writeResult(apply->GetOutput(2), outPrefix + "B1.nii");
	QI::writeResiduals(apply->GetResidOutput(), outPrefix, all_residuals);

	if (verbose) cout << "Finished." << endl;
	return EXIT_SUCCESS;
}
