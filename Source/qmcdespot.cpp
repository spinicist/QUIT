/*
 *  mcdespot_main.cpp
 *
 *  Created by Tobias Wood on 14/02/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <atomic>
#include <getopt.h>
#include <time.h>
#include <fstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "Nifti/Nifti.h"
#include "QUIT/QUIT.h"
#include "Models.h"
#include "Sequence.h"
#include "RegionContraction.h"

using namespace std;
using namespace Eigen;
using namespace QUIT;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
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
	--1, --2, --3     : Use 1, 2 or 3 component sequences (default 3)\n\
	--f0, -f file     : Use f0 Map file (in Hertz)\n\
	--B1, -b file     : B1 Map file (ratio)\n\
	--start, -s n     : Only start processing at slice n.\n\
	--stop, -p n      : Finish at slice n-1\n\
	--scale, -S MEAN  : Normalise signals to mean (default)\n\
	            NONE  : Fit a scaling factor/proton density\n\
	            x     : Fix to x\n\
	--gauss, -g       : Use Gaussian Region Contraction\n\
	--flip, -F        : Data order is phase, then flip-angle (default opposite)\n\
	--tesla, -t 3     : Boundaries suitable for 3T (default)\n\
	            7     : Boundaries suitable for 7T \n\
	            u     : User specified boundaries from stdin\n\
	--sequences, -M s : Use simple sequences (default)\n\
	            f     : Use Finite Pulse Length correction\n\
	--complex, -x     : Fit to complex data\n\
	--contract, -c n  : Read contraction settings from stdin (Will prompt)\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

static shared_ptr<Model> model = make_shared<MCD3>();
static auto tesla = FieldStrength::Three;
static size_t start_slice = 0, stop_slice = numeric_limits<size_t>::max();
static int verbose = false, prompt = true, all_residuals = false,
           fitFinite = false, fitComplex = false, flipData = false,
           gauss = false, samples = 5000, retain = 50, contract = 10,
           voxI = 0, voxJ = 0;
static double expand = 0., scaling = 0.;
static string outPrefix;
static const struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"f0", required_argument, 0, 'f'},
	{"B1", required_argument, 0, 'b'},
	{"start", required_argument, 0, 's'},
	{"stop", required_argument, 0, 'p'},
	{"scale", required_argument, 0, 'S'},
	{"gauss", no_argument, 0, 'g'},
	{"flip", required_argument, 0, 'F'},
	{"tesla", required_argument, 0, 't'},
	{"sequences", no_argument, 0, 'M'},
	{"complex", no_argument, 0, 'x'},
	{"contract", no_argument, 0, 'c'},
	{"resids", no_argument, 0, 'r'},
	{"threads", required_argument, 0, 'T'},
	{"no-prompt", no_argument, 0, 'n'},
	{"1", no_argument, 0, '1'},
	{"2", no_argument, 0, '2'},
	{"3", no_argument, 0, '3'},
	{0, 0, 0, 0}
};
static const char* short_options = "hvm:o:f:b:s:p:S:gt:FT:M:xcrn123i:j:";

//******************************************************************************
#pragma mark Read in all required files and data from cin
//******************************************************************************
Nifti::Header parseInput(SequenceGroup &seq, vector<MultiArray<complex<float>, 4>> &signalVols, Array2d &f0Bandwidth);
Nifti::Header parseInput(SequenceGroup &seq, vector<MultiArray<complex<float>, 4>> &signalVols, Array2d &f0Bandwidth)
{
	Nifti::Header hdr;
	string type, path;
	if (prompt) cout << "Specify next image type (SPGR/SSFP): " << flush;
	f0Bandwidth = Array2d::Zero();
	while (Read(cin, type) && (type != "END") && (type != "")) {
		if (type != "SPGR" && type != "SSFP") {
			throw(std::runtime_error("Unknown signal type: " + type));
		}
		if (prompt) cout << "Enter image path: " << flush;
		Read(cin, path);
		Nifti::File inFile(path);
		if (signalVols.size() == 0) {
			hdr = inFile.header(); // Save header info for later
		} else {
			checkHeaders(hdr, {inFile});
		}
		if (verbose) cout << "Opened: " << inFile.imagePath() << endl;
		Agilent::ProcPar pp; ReadPP(inFile, pp);
		if ((type == "SPGR") && !fitFinite) {
			seq.addSequence(make_shared<SPGRSimple>(prompt, pp));
		} else if ((type == "SPGR" && fitFinite)) {
			seq.addSequence(make_shared<SPGRFinite>(prompt, pp));
		} else if ((type == "SSFP" && !fitFinite)) {
			auto s = make_shared<SSFPSimple>(prompt, pp);
			f0Bandwidth = s->bandwidth();
			seq.addSequence(s);
		} else if ((type == "SSFP" && fitFinite)) {
			auto s = make_shared<SSFPFinite>(prompt, pp);
			f0Bandwidth = s->bandwidth();
			seq.addSequence(s);
		}
		if (seq.sequence(seq.count() - 1)->size() != inFile.dim(4)) {
			throw(std::runtime_error("Number of volumes in file " + inFile.imagePath() + " does not match input."));
		}
		MultiArray<complex<float>, 4> inData(inFile.dims().head(4));
		if (verbose) cout << "Reading data..." << flush;
		inFile.readVolumes(inData.begin(), inData.end());
		signalVols.push_back(inData);
		inFile.close();
		if (verbose) cout << "done." << endl;
		// Print message ready for next loop
		if (prompt) cout << "Specify next image type (SPGR/SSFP, END to finish input): " << flush;
	}
	return hdr;
}

class MCDFunctor : public DenseFunctor<double> {
	public:
		const SequenceBase &m_sequence;
		const ArrayXcd &m_data;
		const bool m_complex, m_debug;
		const shared_ptr<Model> m_model;

		MCDFunctor(SequenceBase &s, const ArrayXcd &d, shared_ptr<Model> m,
		           const bool fitComplex, const bool debug = false) :
			DenseFunctor<double>(m->nParameters(), s.size()),
			m_sequence(s), m_data(d), m_model(m), m_complex(fitComplex), m_debug(debug)
		{
			assert(static_cast<size_t>(m_data.rows()) == values());
		}

		const bool constraint(const VectorXd &params) const {
			return m_model->ValidParameters(params);
		}

		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
			eigen_assert(diffs.size() == values());
			ArrayXcd s = m_sequence.signal(m_model, params);
			if (m_complex) {
				diffs = (s - m_data).abs();
			} else {
				diffs = s.abs() - m_data.abs();
			}
			if (m_debug) {
				cout << endl << __PRETTY_FUNCTION__ << endl;
				cout << "p:     " << params.transpose() << endl;
				cout << "s:     " << s.transpose() << endl;
				cout << "data:  " << m_data.transpose() << endl;
				cout << "diffs: " << diffs.transpose() << endl;
			}
			return 0;
		}
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	try { // To fix uncaught exceptions on Mac
	cout << version << endl << credit_me << endl;
	Eigen::initParallel();
	Nifti::File maskFile, f0File, B1File;
	MultiArray<int8_t, 3> maskVol;
	MultiArray<float, 3> f0Vol, B1Vol;
	ThreadPool threads;

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case '1': model = make_shared<SCD>(); break;
			case '2': model = make_shared<MCD2>(); break;
			case '3': model = make_shared<MCD3>(); break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskVol.resize(maskFile.matrix());
				maskFile.readVolumes(maskVol.begin(), maskVol.end(), 0, 1);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'f':
				cout << "Reading f0 file: " << optarg << endl;
				f0File.open(optarg, Nifti::Mode::Read);
				f0Vol.resize(f0File.matrix());
				f0File.readVolumes(f0Vol.begin(), f0Vol.end(), 0, 1);
				break;
			case 'b':
				cout << "Reading B1 file: " << optarg << endl;
				B1File.open(optarg, Nifti::Mode::Read);
				B1Vol.resize(B1File.matrix());
				B1File.readVolumes(B1Vol.begin(), B1Vol.end(), 0, 1);
				break;
			case 's': start_slice = atoi(optarg); break;
			case 'p': stop_slice = atoi(optarg); break;
			case 'S': {
				string mode(optarg);
				if (mode == "MEAN") {
					if (verbose) cout << "Mean scaling selected." << endl;
					model->setScaling(Model::Scale::ToMean);
				} else if (mode == "FIT") {
					if (verbose) cout << "Fit PD/M0 selected." << endl;
					model->setScaling(Model::Scale::None);
				} else {
					if (verbose) cout << "Scale factor = " << atof(optarg) << endl;
					model->setScaling(Model::Scale::None);
					scaling = atof(optarg);
				}
			} break;
			case 'g':
				gauss = true;
				break;
			case 'F':
				flipData = true;
				break;
			case 'T':
				threads.resize(atoi(optarg));
				break;
			case 't':
				switch (*optarg) {
					case '3': tesla = FieldStrength::Three; break;
					case '7': tesla = FieldStrength::Seven; break;
					case 'u': tesla = FieldStrength::User; break;
					default:
						cout << "Unknown boundaries type " << *optarg << endl;
						return EXIT_FAILURE;
						break;
			} break;
			case 'M':
				switch (*optarg) {
					case 's': fitFinite = false; if (verbose) cout << "Simple sequences selected." << endl; break;
					case 'f': fitFinite = true;  if (verbose) cout << "Finite pulse correction selected." << endl; break;
					default:
						cout << "Unknown sequences type " << *optarg << endl;
						return EXIT_FAILURE;
						break;
			} break;
			case 'x':
				fitComplex = true;
				break;
			case 'c':
				cout << "Enter max number of contractions: " << flush; cin >> contract;
				cout << "Enter number of samples per contraction: " << flush; cin >> samples;
				cout << "Enter number of samples to retain: " << flush; cin >> retain;
				cout << "Enter fraction to expand region by: " << flush; cin >> expand;
				{ string dummy; getline(cin, dummy); } // Eat newlines
				break;
			case 'r': all_residuals = true; break;
			case 'i': voxI = atoi(optarg); break;
			case 'j': voxJ = atoi(optarg); break;
			case 'h':
			case '?': // getopt will print an error message
			default:
				cout << usage << endl;
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 0) {
		cerr << usage << endl << "Incorrect number of arguments." << endl;
		return EXIT_FAILURE;
	}

	//**************************************************************************
	#pragma mark  Read input and set up corresponding SPGR & SSFP lists
	//**************************************************************************
	SequenceGroup sequences;
	Array2d f0Bandwidth;
	// Build a Functor here so we can query number of parameters etc.
	if (verbose) cout << "Using " << model->Name() << " model." << endl;
	vector<MultiArray<complex<float>, 4>> signalVols;
	Nifti::Header hdr = parseInput(sequences, signalVols, f0Bandwidth);
	checkHeaders(hdr, {maskFile, f0File, B1File});
	//**************************************************************************
	#pragma mark Allocate memory and set up boundaries.
	//**************************************************************************
	MultiArray<float, 4> paramsVols(hdr.matrix(), model->nParameters());
	MultiArray<float, 4> ResidsVols(hdr.matrix(), sequences.size());;
	MultiArray<float, 3> ResVol(hdr.matrix());
	
	ArrayXd threshes(model->nParameters()); threshes.setConstant(0.05);
	ArrayXXd bounds = model->Bounds(tesla, 0);
	if (tesla == FieldStrength::User) {
		if (prompt) cout << "Enter parameter pairs (low then high)" << endl;
		for (size_t i = 0; i < model->nParameters() - 1; i++) {
			if (prompt) cout << model->Names()[i] << ": " << flush;
			cin >> bounds(i, 0) >> bounds(i, 1);
		}
	}
	bounds.row(model->nParameters() - 2) = f0Bandwidth;
	ArrayXd weights(sequences.size()); weights.setOnes();
	if (verbose) {
		cout << sequences;
		cout << "Bounds:" << endl <<  bounds.transpose() << endl;
		ofstream boundsFile(outPrefix + "bounds.txt");
		for (size_t p = 0; p < model->nParameters(); p++) {
			boundsFile << model->Names()[p] << "\t" << bounds.row(p) << endl;
		}
		boundsFile.close();
	}
	
	//**************************************************************************
	#pragma mark Do the fitting
	//**************************************************************************
	if (stop_slice > hdr.dim(3))
		stop_slice = hdr.dim(3);
    time_t procStart = time(NULL);
	char theTime[512];
	strftime(theTime, 512, "%H:%M:%S", localtime(&procStart));
	cout << "Started processing at " << theTime << endl;
	for (size_t k = start_slice; k < stop_slice; k++) {
		if (verbose) cout << "Processing slice " << k << "..." << flush;
		atomic<int> voxCount{0};
		clock_t loopStart = clock();

		function<void (const size_t, const size_t)>
		processVox = [&] (const size_t i, const size_t j) {
			if (!maskFile || maskVol[{i,j,k}]) {
				ArrayXcd signal = sequences.loadSignals(signalVols, i, j, k, flipData);
				double B1 = B1File ? B1Vol[{i,j,k}] : 1.;
				ArrayXXd localBounds = bounds;
				if (f0File) {
					localBounds.row(model->nParameters() - 2).setConstant(f0Vol[{i,j,k}]);
				}
				localBounds.row(model->nParameters() - 1).setConstant(B1);
				if (model->scaling() == Model::Scale::None) {
					if (scaling == 0.) {
						localBounds(0, 0) = 0.;
						localBounds(0, 1) = signal.abs().maxCoeff() * 25;
					} else {
						localBounds.row(0).setConstant(scaling);
					}
				}
				MCDFunctor func(sequences, signal, model, fitComplex, false);
				RegionContraction<MCDFunctor> rc(func, localBounds, weights, threshes,
													samples, retain, contract, expand, gauss, (voxI != 0));
				ArrayXd params(model->nParameters());
				rc.optimise(params); // Add the voxel number to the time to get a decent random seed
				paramsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = params.cast<float>();
				ResidsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = rc.residuals().cast<float>();
				ResVol[{i,j,k}] = static_cast<float>(rc.SoS());
				if ((rc.status() == RCStatus::Converged) || (rc.status() == RCStatus::IterationLimit)) {
					voxCount++;
				}
			}
		};
		if (voxI == 0) {
			threads.for_loop2(processVox, hdr.dim(1), hdr.dim(2));
			if (threads.interrupted())
				break;
		} else {
			processVox(voxI, voxJ);
			return EXIT_SUCCESS;
		}
		if (verbose) {
			clock_t loopEnd = clock();
			if (voxCount > 0)
				cout << voxCount << " unmasked voxels, CPU time per voxel was "
				          << ((loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC)) << " s, ";
			cout << "finished." << endl;
		}
		if (threads.interrupted())
			break;
	}
	time_t procEnd = time(NULL);
	strftime(theTime, 512, "%H:%M:%S", localtime(&procEnd));
	cout << "Finished processing at " << theTime << ". Run-time was " 
		 << difftime(procEnd, procStart) << " s." << endl;
	// Residuals can only be written here if we want them to go in a 4D gzipped file
	outPrefix = outPrefix + model->Name() + "_";
	hdr.setDim(4, 1);
	hdr.setDatatype(Nifti::DataType::FLOAT32);
	hdr.description = version;
	hdr.intent = Nifti::Intent::Estimate;
	size_t start = (model->scaling() == Model::Scale::None) ? 0 : 1; // Skip PD
	size_t end = model->nParameters() - 1; // Skip B1 for now
	for (size_t p = start; p < end; p++) {
		hdr.intent_name = model->Names().at(p);
		Nifti::File file(hdr, outPrefix + model->Names().at(p) + "" + OutExt());
		auto param = paramsVols.slice<3>({0,0,0,p},{-1,-1,-1,0});
		file.writeVolumes(param.begin(), param.end());
		file.close();
	}
	hdr.intent_name = "Fractional Residual";
	Nifti::File SoS(hdr, outPrefix + "residual" + OutExt());
	SoS.writeVolumes(ResVol.begin(), ResVol.end());
	SoS.close();
	if (all_residuals) {
		hdr.setDim(4, static_cast<int>(sequences.size()));
		hdr.intent_name = "Residuals";
		Nifti::File res(hdr, outPrefix + "residuals" + OutExt());
		res.writeVolumes(ResidsVols.begin(), ResidsVols.end());
		res.close();
	}
	cout << "Finished writing data." << endl;
	
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

