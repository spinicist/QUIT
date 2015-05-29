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

#include "Sequence.h"
#include "Util.h"

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
	--non-ge, -N      : Use a generic MP-RAGE sequence, not GE IR-SPGR\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--mask, -m file   : Mask input with specified file\n\
	--thresh, -t n    : Threshold maps at PD < n\n\
	--clamp, -c n     : Clamp T1 between 0 and n\n\
	--start, -s N     : Start processing from slice N\n\
	--stop, -p  N     : Stop processing at slice N\n\
	--its, -i N       : Max iterations for NLLS (default 4)\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

static bool verbose = false, prompt = true, GE = true;
static size_t nIterations = 4;
static string outPrefix;
static double thresh = -numeric_limits<double>::infinity();
static double clamp_lo = -numeric_limits<double>::infinity(), clamp_hi = numeric_limits<double>::infinity();
static size_t start_slice = 0, stop_slice = numeric_limits<size_t>::max();
static const struct option long_opts[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"non-ge", no_argument, 0, 'N'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"thresh", required_argument, 0, 't'},
	{"clamp", required_argument, 0, 'c'},
	{"start", required_argument, 0, 's'},
	{"stop", required_argument, 0, 'p'},
	{"its", required_argument, 0, 'i'},
	{"threads", required_argument, 0, 'T'},
	{0, 0, 0, 0}
};
static const char *short_opts = "hvnNm:o:t:c:s:p:i:T:";

// HIFI Algorithm - includes optimising B1
class HIFIFunctor : public DenseFunctor<double> {
	protected:
		const SequenceBase &m_sequence;
		const ArrayXd m_data;
		const shared_ptr<SCD> m_model = make_shared<SCD>();

	public:
		HIFIFunctor(SequenceBase &cs, const ArrayXd &data, const bool debug) :
			DenseFunctor<double>(3, cs.size()),
			m_sequence(cs), m_data(data), m_debug(debug)
		{
			assert(static_cast<size_t>(m_data.rows()) == values());
		}

		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
			eigen_assert(diffs.size() == values());
			VectorXd fullpar(5); // Add T2 and f0
			fullpar << params(0), params(1), 0, 0, params(2);
			ArrayXcd s = m_sequence.signal(m_model, fullpar);
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
			outputs << data.maxCoeff() * 10., 1., 1.; // Initial guess
			lm.lmder1(outputs);
			VectorXd pfull(5); pfull << outputs[0], outputs[1], 0, 0, outputs[2]; // Build full parameter vector
			ArrayXd theory = combined.signal(model, pfull).abs();
			resids = (data.array() - theory);
		}
};


//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	cout << version << endl << "Improved formulas thanks to Michael Thrippleton." << endl;
	Eigen::initParallel();
	Nifti::File maskFile, spgrFile, irFile;
	MultiArray<int8_t, 3> maskVol;
	ThreadPool threads;
	shared_ptr<SCD> model = make_shared<SCD>();

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'N': GE = false; break;
			case 'm':
				if (verbose) cout << "Opening mask file: " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskVol.resize(maskFile.matrix());
				maskFile.readVolumes(maskVol.begin(), maskVol.end(), 0, 1);
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
			case 's': start_slice = atoi(optarg); break;
			case 'p': stop_slice = atoi(optarg); break;
			case 'i':
				nIterations = atoi(optarg);
				break;
			case 'T':
				threads.resize(atoi(optarg));
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
	
	SequenceGroup combined;
	if (verbose) cout << "Opening SPGR file: " << argv[optind] << endl;
	spgrFile.open(argv[optind++], Nifti::Mode::Read);
	Agilent::ProcPar pp; ReadPP(spgrFile, pp);
	shared_ptr<SteadyState> spgrSequence = make_shared<SPGRSimple>(prompt, pp);
	combined.addSequence(spgrSequence);

	if (verbose) cout << "Opening IR-SPGR file: " << argv[optind] << endl;
	irFile.open(argv[optind], Nifti::Mode::Read);
	shared_ptr<SteadyState> irspgr;
	if (GE) {
		irspgr = make_shared<IRSPGR>(prompt);
	} else {
		irspgr = make_shared<MPRAGE>(prompt);
	}
	combined.addSequence(irspgr);

	if (spgrSequence->size() != spgrFile.header().dim(4)) {
		throw(std::runtime_error("Specified number of flip-angles does not match number of volumes in file: " + spgrFile.imagePath()));
	}
	if (irspgr->size() != irFile.header().dim(4)) {
		throw(std::runtime_error("Specified number of TIs does not match number of volumes in file: " + irFile.imagePath()));
	}
	checkHeaders(spgrFile.header(),{irFile, maskFile});
	if (verbose) {
		cout << combined << endl;
		cout << "Reading image data..." << flush;
	}
	const auto dims = spgrFile.matrix();
	MultiArray<float, 4> SPGR_Vols(dims, spgrFile.dim(4));
	MultiArray<float, 4> IR_Vols(dims, irFile.dim(4));
	spgrFile.readVolumes(SPGR_Vols.begin(), SPGR_Vols.end());
	irFile.readVolumes(IR_Vols.begin(), IR_Vols.end());
	spgrFile.close();
	irFile.close();
	if (verbose) cout << "done" << endl;

	MultiArray<float, 3> PD_Vol(dims), T1_Vol(dims), B1_Vol(dims), res_Vol(dims);
	if (stop_slice > dims[2])
		stop_slice = dims[2];
	time_t startTime;
	if (verbose)
		startTime = printStartTime();
	clock_t startClock = clock();
	int voxCount = 0;
	for (size_t k = start_slice; k < stop_slice; k++) {
		clock_t loopStart = clock();
		// Read in data
		if (verbose) cout << "Starting slice " << k << "..." << flush;
		atomic<int> sliceCount{0};
		function<void (const size_t&, const size_t&)> processVox = [&] (const size_t &j, const size_t &i) {
			double T1 = 0., PD = 0., B1 = 1., res = 0.; // Assume B1 field is uniform for classic DESPOT
			const MultiArray<float, 3>::Index idx{i,j,k};
			if (!maskFile || (maskVol[idx] > 0.)) {
				sliceCount++;
				ArrayXd spgrSig = SPGR_Vols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray().abs().cast<double>();

				// Get a first guess with DESPOT1
				VectorXd Y = spgrSig / spgrSequence->flip().sin();
				MatrixXd X(Y.rows(), 2);
				X.col(0) = spgrSig / spgrSequence->flip().tan();
				X.col(1).setOnes();
				VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
				T1 = -spgrSequence->TR() / log(b[0]);
				PD = b[1] / (1. - b[0]);

				ArrayXd irSig = IR_Vols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray().abs().cast<double>();
				ArrayXd combinedSig(combined.size());
				combinedSig.head(spgrSig.size()) = spgrSig;
				combinedSig.tail(irSig.size()) = irSig;

				// Now add in IR-SPGR data and do a Lev-Mar fit
				HIFIFunctor f(combined, combinedSig, false);
				NumericalDiff<HIFIFunctor> nDiff(f);
				LevenbergMarquardt<NumericalDiff<HIFIFunctor>> lm(nDiff);
				lm.setMaxfev(nIterations);
				VectorXd p(3);
				p << PD, T1, B1; // Don't need T2 of f0 for this (yet)
				lm.lmder1(p);
				PD = p(0); T1 = p(1); B1 = p(2);
				if (PD < thresh) {
					PD = 0.;
					T1 = 0.;
					B1 = 0.;
				}
				T1 = clamp(T1, clamp_lo, clamp_hi);
				VectorXd pfull(5); pfull << PD, T1, 0, 0, B1;
				ArrayXd theory = combined.signal(model, pfull).abs();
				ArrayXd resids = (combinedSig - theory);
				res = sqrt(resids.square().sum() / resids.rows()) / PD;
			}
			PD_Vol[idx] = PD;
			T1_Vol[idx] = T1;
			B1_Vol[idx] = B1;
			res_Vol[idx] = res;
		};
		threads.for_loop2(processVox, dims(1), dims(0));
		if (verbose) printLoopTime(loopStart, sliceCount);
		voxCount += sliceCount;
		if (threads.interrupted())
			break;

	}
	if (verbose) {
		printElapsedTime(startTime);
		printElapsedClock(startClock, voxCount);
		cout << "Writing results." << endl;
	}
	//**************************************************************************
	#pragma mark Write out data
	//**************************************************************************
	Nifti::Header outHdr = spgrFile.header();
	outPrefix = outPrefix + "HIFI_";
	outHdr.description = version;
	outHdr.setDim(4, 1);
	outHdr.setDatatype(Nifti::DataType::FLOAT32);
	outHdr.intent = Nifti::Intent::Estimate;
	outHdr.intent_name = "T1 (seconds)";
	Nifti::File outFile(outHdr, outPrefix + "T1" + OutExt());
	outFile.writeVolumes(T1_Vol.begin(), T1_Vol.end());
	outFile.close();
	outHdr.intent_name = "PD (au)";
	outFile.setHeader(outHdr);
	outFile.open(outPrefix + "PD" + OutExt(), Nifti::Mode::Write);
	outFile.writeVolumes(PD_Vol.begin(), PD_Vol.end());
	outFile.close();
	outHdr.intent_name = "B1 Field Ratio";
	outFile.setHeader(outHdr);
	outFile.open(outPrefix + "B1" + OutExt(), Nifti::Mode::Write);
	outFile.writeVolumes(B1_Vol.begin(), B1_Vol.end());
	outFile.close();
	outHdr.intent_name = "Fractional Residual";
	outFile.setHeader(outHdr);
	outFile.open(outPrefix + "residual" + OutExt(), Nifti::Mode::Write);
	outFile.writeVolumes(res_Vol.begin(), res_Vol.end());
	outFile.close();
	if (verbose) cout << "Finished." << endl;
	return EXIT_SUCCESS;
}
