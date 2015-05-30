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
using namespace QUITK;

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
		size_t m_iterations = 4;
	public:
		void setIterations(size_t n) { m_iterations = n; }

		size_t numConsts() const override { return 0; }
		size_t numOutputs() const override { return 3; }

		virtual VectorXd defaultConsts() {
			// No constants for HIFI
			VectorXd def = VectorXd::Zero(0);
			return def;
		}

		virtual void apply(const shared_ptr<SequenceBase> sequence,
		                   const VectorXd &data,
		                   const VectorXd &inputs,
		                   VectorXd &outputs,
		                   ArrayXd &resids) const override
		{
			HIFIFunctor f(sequence, data);
			NumericalDiff<HIFIFunctor> nDiff(f);
			LevenbergMarquardt<NumericalDiff<HIFIFunctor>> lm(nDiff);
			lm.setMaxfev(m_iterations * (sequence->size() + 1));
			// PD, T1, B1
			outputs << data.array().abs().maxCoeff() * 50., 1., 1.; // Initial guess
			lm.lmder1(outputs);
			VectorXd pfull(5); pfull << outputs[0], outputs[1], 0, 0, outputs[2]; // Build full parameter vector
			auto model = make_shared<SCD>();
			ArrayXd theory = sequence->signal(model, pfull).abs();
			resids = (data.array() - theory);
		}
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();

	typedef itk::Image<float, 3> FloatImage;
	typedef itk::VectorImage<float, 3> FloatVectorImage;

	typedef itk::ImageFileReader<FloatImage> Reader;
	typedef itk::ImageFileReader<itk::Image<float, 4>> Reader4D;
	Reader::Pointer mask = ITK_NULLPTR;
	Reader::Pointer B1   = ITK_NULLPTR;

	auto hifi = make_shared<HIFIAlgo>();
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'M': IR = false; break;
			case 'm':
				if (verbose) cout << "Opening mask file: " << optarg << endl;
				mask = Reader::New();
				mask->SetFileName(optarg);
				break;
			case 'o':
				outPrefix = optarg;
				if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 't': thresh = atof(optarg); break;
			case 'c':
				clamp_lo = 0;
				clamp_hi = atof(optarg);
				break;
			case 'i':
				hifi->setIterations(atoi(optarg));
				break;
			case 'r': all_residuals = true; break;
			case 'T':
				itk::MultiThreader::SetGlobalDefaultNumberOfThreads(atoi(optarg));
				break;
			case 'h':
			case '?': // getopt will print an error message
				cout << usage << endl;
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 2) {
		cerr << "Incorrect number of arguments." << endl;
		cout << usage << endl;
		return EXIT_FAILURE;
	}
	
	if (verbose) cout << "Opening SPGR file: " << argv[optind] << endl;
	auto spgrFile = Reader4D::New();
	auto spgrImg = ImageToVectorFilter<float>::New();
	spgrFile->SetFileName(argv[optind++]);
	spgrImg->SetInput(spgrFile->GetOutput());
	auto spgrSequence = make_shared<SPGRSimple>(prompt);

	if (verbose) cout << "Opening IR-SPGR file: " << argv[optind] << endl;
	auto irFile = Reader4D::New();
	auto irImg = ImageToVectorFilter<float>::New();
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
	auto apply = itk::ApplyAlgorithmFilter<float, HIFIAlgo>::New();
	apply->SetSequence(combined);
	apply->SetAlgorithm(hifi);
	apply->Setup();
	apply->SetDataInput(0, spgrImg->GetOutput());
	apply->SetDataInput(1, irImg->GetOutput());
	if (mask)
		apply->SetMask(mask->GetOutput());
	if (verbose) cout << "Processing..." << endl;
	apply->Update();
	if (verbose) cout << "Writing results." << endl;
	outPrefix = outPrefix + "HIFI_";

	QUITK::writeResult(apply->GetOutput(0), outPrefix + "PD.nii");
	QUITK::writeResult(apply->GetOutput(1), outPrefix + "T1.nii");
	QUITK::writeResult(apply->GetOutput(2), outPrefix + "B1.nii");
	QUITK::writeResiduals(apply->GetResidOutput(), outPrefix, all_residuals);

	if (verbose) cout << "Finished." << endl;
	return EXIT_SUCCESS;
}
