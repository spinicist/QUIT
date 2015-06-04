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

#include "Util.h"
#include "Filters/ApplyAlgorithmFilter.h"
#include "Filters/ReorderVectorFilter.h"
#include "Model.h"
#include "Sequence.h"
#include "RegionContraction.h"

using namespace std;
using namespace Eigen;

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
	--contract, -c n  : Read contraction settings from stdin (Will prompt)\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

static auto tesla = FieldStrength::Three;
static size_t start_slice = 0, stop_slice = numeric_limits<size_t>::max();
static int verbose = false, prompt = true, all_residuals = false,
           fitFinite = false, flipData = false;
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
	/*{"complex", no_argument, 0, 'x'},*/
	{"contract", no_argument, 0, 'c'},
	{"resids", no_argument, 0, 'r'},
	{"threads", required_argument, 0, 'T'},
	{"no-prompt", no_argument, 0, 'n'},
	{"1", no_argument, 0, '1'},
	{"2", no_argument, 0, '2'},
	{"3", no_argument, 0, '3'},
	{0, 0, 0, 0}
};
static const char* short_options = "hvm:o:f:b:s:p:S:gt:FT:M:crn123i:j:";

/*
 * Read in all required files and data from cin
 */
void parseInput(shared_ptr<SequenceGroup> seq, vector<typename QUITK::ReadFloatTimeseries::Pointer> &files, vector<typename QUITK::FloatTimeseriesToVector::Pointer> &data, vector<typename itk::ReorderVectorFilter<QUITK::FloatVectorImage>::Pointer> &order, Array2d &f0Bandwidth, bool flip);
void parseInput(shared_ptr<SequenceGroup> seq,
                vector<typename QUITK::ReadFloatTimeseries::Pointer> &files,
                vector<typename QUITK::FloatTimeseriesToVector::Pointer> &data,
                vector<typename itk::ReorderVectorFilter<QUITK::FloatVectorImage>::Pointer> &order,
                Array2d &f0Bandwidth, bool flip)
{
	string type, path;
	if (prompt) cout << "Specify next image type (SPGR/SSFP): " << flush;
	f0Bandwidth = Array2d::Zero();
	while (QUITK::Read(cin, type) && (type != "END") && (type != "")) {
		if (type != "SPGR" && type != "SSFP") {
			throw(std::runtime_error("Unknown signal type: " + type));
		}
		if (prompt) cout << "Enter image path: " << flush;
		QUITK::Read(cin, path);
		files.push_back(QUITK::ReadFloatTimeseries::New());
		files.back()->SetFileName(path);
		data.push_back(QUITK::FloatTimeseriesToVector::New());
		data.back()->SetInput(files.back()->GetOutput());
		order.push_back(itk::ReorderVectorFilter<QUITK::FloatVectorImage>::New());
		order.back()->SetInput(data.back()->GetOutput());
		if (verbose) cout << "Opened: " << path << endl;
		if ((type == "SPGR") && !fitFinite) {
			seq->addSequence(make_shared<SPGRSimple>(prompt));
		} else if ((type == "SPGR" && fitFinite)) {
			seq->addSequence(make_shared<SPGRFinite>(prompt));
		} else if ((type == "SSFP" && !fitFinite)) {
			auto s = make_shared<SSFPSimple>(prompt);
			f0Bandwidth = s->bandwidth();
			seq->addSequence(s);
			if (flip)
				order.back()->SetStride(s->phases());
		} else if ((type == "SSFP" && fitFinite)) {
			auto s = make_shared<SSFPFinite>(prompt);
			f0Bandwidth = s->bandwidth();
			seq->addSequence(s);
			if (flip)
				order.back()->SetStride(s->phases());
		}
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
			diffs = s.abs() - m_data;
			return 0;
		}
};

class MCDAlgo : public Algorithm<double> {
	private:
		size_t m_samples = 5000, m_retain = 50, m_contractions = 10;
		bool m_gauss = false;
		ArrayXXd m_bounds;
		shared_ptr<Model> m_model = nullptr;

	public:
		void setModel(shared_ptr<Model> &m) { m_model = m; }
		void setSamples(size_t s) { m_samples = s; }
		void setRetain(size_t r) { m_retain = r; }
		void setRCPars(size_t c, size_t s, size_t r) { m_contractions = c; m_samples = s; m_retain = r; }
		void setScaling(Model::Scale s) { m_model->setScaling(s); }
		void setBounds(ArrayXXd b) { m_bounds = b; }
		void setGauss(bool g) { m_gauss = g; }

		size_t numConsts() const override { return 2; }
		size_t numOutputs() const override { return m_model->nParameters(); }

		virtual VectorXd defaultConsts() {
			// f0, B1
			VectorXd def = VectorXd::Ones(2);
			def[0] = NAN;
			return def;
		}

		virtual void apply(const shared_ptr<SequenceBase> sequence,
		                   const VectorXd &data,
		                   const VectorXd &inputs,
		                   VectorXd &outputs,
		                   ArrayXd &resids) const override
		{
			ArrayXd thresh(m_model->nParameters()); thresh.setConstant(0.05);
			ArrayXd weights = ArrayXd::Ones(sequence->size());
			double f0 = inputs[0];
			double B1 = inputs[1];
			ArrayXXd localBounds = m_bounds;
			if (!std::isnan(f0))
				localBounds.row(m_model->nParameters() - 2).setConstant(f0);
			if (m_model->scaling() == Model::Scale::None) {
				if (scaling == 0.) {
					localBounds(0, 0) = 0.;
					localBounds(0, 1) = data.array().maxCoeff() * 25;
				} else {
					localBounds.row(0).setConstant(scaling);
				}
			}
			localBounds.row(m_model->nParameters() - 1).setConstant(B1);
			//cout << "T1 " << T1 << " B1 " << B1 << " inputs " << inputs.transpose() << endl;
			//cout << bounds << endl;
			//cout << "data " << data.transpose() << endl;
			//cout << "f0 " << f0 << " B1 " << B1 << endl;
			//cout << "lb" << endl << localBounds.transpose() << endl;
			//cout << "data " << data.transpose() << endl;
			MCDFunctor func(m_model, sequence, data);
			//cout << "data " << data.transpose() << endl;
			RegionContraction<MCDFunctor> rc(func, localBounds, weights, thresh,
			                                m_samples, m_retain, m_contractions, 0., true, false);
			//cout << "data " << data.transpose() << endl;
			rc.optimise(outputs);
			//cout << "outputs " << outputs.transpose() << endl;
			resids = rc.residuals();
		}
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();
	QUITK::ReadFloatImage::Pointer mask, B1, f0 = ITK_NULLPTR;
	shared_ptr<MCDAlgo> mcd = make_shared<MCDAlgo>();
	shared_ptr<Model> model;
	// Deal with these options in first pass to ensure the correct model is selected
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case '1': model = make_shared<SCD>(); break;
			case '2': model = make_shared<MCD2>(); break;
			case '3': model = make_shared<MCD3>(); break;
			default:
				break;
		}
	}
	// Now reset and do a second pass
	optind = 1;
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'm':
				if (verbose) cout << "Reading mask file " << optarg << endl;
				mask = QUITK::ReadFloatImage::New();
				mask->SetFileName(optarg);
				break;
			case 'o':
				outPrefix = optarg;
				if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'f':
				if (verbose) cout << "Reading f0 file: " << optarg << endl;
				f0 = QUITK::ReadFloatImage::New();
				f0->SetFileName(optarg);
				break;
			case 'b':
				if (verbose) cout << "Reading B1 file: " << optarg << endl;
				B1 = QUITK::ReadFloatImage::New();
				B1->SetFileName(optarg);
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
			case 'g': mcd->setGauss(true); break;
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
			case 'M':
				switch (*optarg) {
					case 's': fitFinite = false; if (verbose) cout << "Simple sequences selected." << endl; break;
					case 'f': fitFinite = true;  if (verbose) cout << "Finite pulse correction selected." << endl; break;
					default:
						cout << "Unknown sequences type " << *optarg << endl;
						return EXIT_FAILURE;
						break;
			} break;
			case 'c': {
				if (prompt) cout << "Enter max contractions/samples per contraction/retained samples/expand fraction: " << flush;
				ArrayXi in = ArrayXi::Zero(3);
				QUITK::ReadEigen(cin, in);
				mcd->setRCPars(in[0], in[1], in[2]);
			} break;
			case 'r': all_residuals = true; break;
			case 'h':
			case '?': // getopt will print an error message
				cout << usage << endl;
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 0) {
		cerr << usage << endl << "Incorrect number of arguments." << endl;
		return EXIT_FAILURE;
	}

	shared_ptr<SequenceGroup> sequences = make_shared<SequenceGroup>();
	Array2d f0Bandwidth;
	// Build a Functor here so we can query number of parameters etc.
	if (verbose) cout << "Using " << model->Name() << " model." << endl;
	vector<QUITK::ReadFloatTimeseries::Pointer> inFiles;
	vector<QUITK::FloatTimeseriesToVector::Pointer> inData;
	vector<itk::ReorderVectorFilter<QUITK::FloatVectorImage>::Pointer> inOrder;
	parseInput(sequences, inFiles, inData, inOrder, f0Bandwidth, flipData);

	ArrayXXd bounds = model->Bounds(tesla, 0);
	if (tesla == FieldStrength::User) {
		if (prompt) cout << "Enter parameter pairs (low then high)" << endl;
		for (size_t i = 0; i < model->nParameters() - 1; i++) {
			if (prompt) cout << model->Names()[i] << ": " << flush;
			cin >> bounds(i, 0) >> bounds(i, 1);
		}
	}
	bounds.row(model->nParameters() - 2) = f0Bandwidth;
	mcd->setModel(model);
	mcd->setBounds(bounds);
	if (verbose) {
		cout << *sequences;
		cout << "Bounds:" << endl <<  bounds.transpose() << endl;
		ofstream boundsFile(outPrefix + "bounds.txt");
		for (size_t p = 0; p < model->nParameters(); p++) {
			boundsFile << model->Names()[p] << "\t" << bounds.row(p) << endl;
		}
		boundsFile.close();
	}
	
	auto apply = itk::ApplyAlgorithmFilter<float, MCDAlgo>::New();
	apply->SetSequence(sequences);
	apply->SetAlgorithm(mcd);
	apply->Setup();
	for (int i = 0; i < inOrder.size(); i++) {
		apply->SetDataInput(i, inOrder.at(i)->GetOutput());
	}
	if (f0)
		apply->SetConstInput(0, f0->GetOutput());
	if (B1)
		apply->SetConstInput(1, B1->GetOutput());
	if (mask)
		apply->SetMask(mask->GetOutput());

	time_t startTime;
	if (verbose) startTime = QUITK::printStartTime();
	apply->Update();
	QUITK::printElapsedTime(startTime);

	outPrefix = outPrefix + model->Name() + "_";
	for (int i = 0; i < model->nParameters(); i++) {
		QUITK::writeResult(apply->GetOutput(i), outPrefix + model->Names()[i] + QUITK::OutExt());
	}
	QUITK::writeResiduals(apply->GetResidOutput(), outPrefix, all_residuals);
	return EXIT_SUCCESS;
}

