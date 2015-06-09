/*
 *  despot2fm.cpp
 *
 *  Created by Tobias Wood on 2015/06/03.
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

#include "Util.h"
#include "Filters/ApplyAlgorithmFilter.h"
#include "Model.h"
#include "Sequence.h"
#include "RegionContraction.h"


using namespace std;
using namespace Eigen;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: despot2-fm [options] T1_map ssfp_file\n\
\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print slice processing times\n\
	--no-prompt, -n   : Suppress input prompts\n\
	--mask, -m file   : Mask input with specified file\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--f0, -f file     : Use f0 Map file (in Hertz)\n\
	--B1, -b file     : B1 Map file (ratio)\n\
	--start, -s N     : Start processing from slice N\n\
	--stop, -p  N     : Stop processing at slice N\n\
	--scale, -S 0     : Normalise signals to mean (default)\n\
	            1     : Fit a scaling factor/proton density\n\
	--flip, -F        : Data order is phase, then flip-angle (default opposite)\n\
	--sequences, -M s : Use simple sequences (default)\n\
	            f     : Use finite pulse length correction\n\
	--contract, -c n  : Read contraction settings from stdin (Will prompt)\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};
/* --complex, -x     : Fit to complex data\n\ */

static auto tesla = FieldStrength::Three;
static int start_slice = 0, stop_slice = 0;
static int verbose = false, prompt = true, all_residuals = false,
           fitFinite = false, fitComplex = false, flipData = false,
           samples = 2000, retain = 20, contract = 10,
           seed = -1;
static double expand = 0.;
static string outPrefix;
static struct option long_opts[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"f0", required_argument, 0, 'f'},
	{"B1", required_argument, 0, 'b'},
	{"start", required_argument, 0, 's'},
	{"stop", required_argument, 0, 'p'},
	{"scale", required_argument, 0, 'S'},
	{"flip", required_argument, 0, 'F'},
	{"threads", required_argument, 0, 'T'},
	{"sequences", no_argument, 0, 'M'},
	{"contract", no_argument, 0, 'c'},
	{"resids", no_argument, 0, 'r'},
	{0, 0, 0, 0}
};
static const char* short_opts = "hvnm:o:f:b:s:p:S:FT:M:crd:";

class FMFunctor : public DenseFunctor<double> {
	public:
		const shared_ptr<SequenceBase> m_sequence;
		shared_ptr<SCD> m_model;
		ArrayXd m_data;
		const double m_T1, m_B1;

		FMFunctor(const shared_ptr<SCD> m, const double T1, const shared_ptr<SequenceBase> s, const ArrayXd &d, const double B1) :
			DenseFunctor<double>(3, s->size()),
			m_model(m), m_sequence(s), m_data(d),
			m_T1(T1), m_B1(B1)
		{
			assert(static_cast<size_t>(m_data.rows()) == values());
		}

		const bool constraint(const VectorXd &params) const {
			Array4d fullparams;
			fullparams << params(0), m_T1, params(1), params(2);
			return m_model->ValidParameters(fullparams);
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

class FMAlgo : public Algorithm<double> {
	private:
		size_t m_samples = 2000, m_retain = 20, m_contractions = 10;
		Array2d m_f0Bounds = Array2d::Zero();
		const shared_ptr<SCD> m_model = make_shared<SCD>();
		shared_ptr<SteadyState> m_sequence;

	public:
		void setSequence(shared_ptr<SteadyState> s) { m_sequence = s; }
		void setSamples(size_t s) { m_samples = s; }
		void setRetain(size_t r) { m_retain = r; }
		void setContractions(size_t c) { m_contractions = c; }
		void setScaling(Model::Scale s) { m_model->setScaling(s); }
		void setf0Bounds(Array2d b) { m_f0Bounds = b; }

		size_t numInputs() const override  { return m_sequence->count(); }
		size_t numConsts() const override  { return 2; }
		size_t numOutputs() const override { return 3; }
		size_t dataSize() const override   { return m_sequence->size(); }

		virtual VectorXd defaultConsts() {
			// T1 & B1
			VectorXd def = VectorXd::Ones(2);
			return def;
		}

		virtual void apply(const VectorXd &data, const VectorXd &inputs,
		                   VectorXd &outputs, ArrayXd &resids) const override
		{
			ArrayXd thresh(3); thresh.setConstant(0.05);
			ArrayXd weights(m_sequence->size()); weights.setOnes();
			ArrayXXd bounds = ArrayXXd::Zero(3, 2);
			double T1 = inputs[0];
			if (isfinite(T1) && (T1 > 0.001)) {
				double B1 = inputs[1];
				if (m_model->scaling() == Model::Scale::None) {
					bounds(0, 0) = 0.;
					bounds(0, 1) = data.array().abs().maxCoeff() * 25;
				} else {
					bounds.row(0).setConstant(1.);
				}
				bounds(1,0) = 0.001;
				bounds(1,1) = T1;
				bounds.row(2) = m_f0Bounds;
				//cout << "T1 " << T1 << " B1 " << B1 << " inputs " << inputs.transpose() << endl;
				//cout << bounds << endl;
				FMFunctor func(m_model, T1, m_sequence, data, B1);
				RegionContraction<FMFunctor> rc(func, bounds, weights, thresh,
												m_samples, m_retain, m_contractions, expand, false, false);
				rc.optimise(outputs);
				resids = rc.residuals();
			} else {
				// No point in processing -ve T1
				outputs.setZero();
				resids.setZero();
			}
		}
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();
	QI::ReadImageF::Pointer mask, B1, f0 = ITK_NULLPTR;
	shared_ptr<FMAlgo> fm = make_shared<FMAlgo>();
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
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
			case 'S':
				switch (atoi(optarg)) {
					case 0 : fm->setScaling(Model::Scale::ToMean); break;
					case 1 : fm->setScaling(Model::Scale::None); break;
					default:
						cout << "Invalid scaling mode: " + to_string(atoi(optarg)) << endl;
						return EXIT_FAILURE;
						break;
				} break;
			case 'F': flipData = true; break;
			case 'T': itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg)); break;
			case 'd': seed = atoi(optarg); break;
			case 'M':
				switch (*optarg) {
					case 's': fitFinite = false; cout << "Simple sequences selected." << endl; break;
					case 'f': fitFinite = true; cout << "Finite pulse correction selected." << endl; break;
					default:
						cout << "Unknown sequences type " << *optarg << endl;
						return EXIT_FAILURE;
						break;
				}
				break;
			case 'c':
				cout << "Enter max number of contractions: " << flush; cin >> contract;
				cout << "Enter number of samples per contraction: " << flush; cin >> samples;
				cout << "Enter number of samples to retain: " << flush; cin >> retain;
				cout << "Enter fraction to expand region by: " << flush; cin >> expand;
				break;
			case 'r': all_residuals = true; break;
			case '?': // getopt will print an error message
			case 'h':
			default:
				cout << usage << endl;
				return EXIT_SUCCESS;
				break;
		}
	}
	if ((argc - optind) != 2) {
		cout << "Wrong number of arguments. Need a T1 map and one SSFP file." << endl;
		return EXIT_FAILURE;
	}
	shared_ptr<SSFPSimple> ssfpSequence;
	if (fitFinite) {
		ssfpSequence = make_shared<SSFPFinite>(prompt);
	} else {
		ssfpSequence = make_shared<SSFPSimple>(prompt);
	}
	if (verbose) cout << *ssfpSequence << endl;

	if (verbose) cout << "Reading T1 Map from: " << argv[optind] << endl;
	auto T1File = QI::ReadImageF::New();
	auto T1Slice = QI::ImageROIF::New(), B1Slice = QI::ImageROIF::New(),
	     maskSlice = QI::ImageROIF::New();

	T1File->SetFileName(argv[optind++]);
	T1File->Update(); // Need to have some image dimensions to do this
	QI::ImageF::RegionType region = T1File->GetOutput()->GetLargestPossibleRegion();
	region.GetModifiableIndex()[2] = start_slice;
	if (stop_slice != 0)
		region.GetModifiableSize()[2] = stop_slice - start_slice;
	else
		region.GetModifiableSize()[2] = region.GetSize()[2] - start_slice;
	T1Slice->SetInput(T1File->GetOutput());
	T1Slice->SetRegionOfInterest(region);
	if (verbose) cout << "Processing region: " << endl << region << endl;
	if (verbose) cout << "Opening SSFP file: " << argv[optind] << endl;
	auto ssfpFile = QI::ReadTimeseriesF::New();
	auto ssfpData = QI::TimeseriesToVectorF::New();
	auto ssfpFlip = QI::ReorderF::New();
	auto ssfpSlice = QI::VectorImageROIF::New();
	ssfpFile->SetFileName(argv[optind++]);
	ssfpData->SetInput(ssfpFile->GetOutput());
	ssfpFlip->SetInput(ssfpData->GetOutput());
	if (flipData) {
		ssfpFlip->SetStride(ssfpSequence->phases());
	}
	ssfpSlice->SetInput(ssfpFlip->GetOutput());
	ssfpSlice->SetRegionOfInterest(region);

	auto apply = itk::ApplyAlgorithmFilter<QI::VectorImageF, FMAlgo>::New();
	fm->setSequence(ssfpSequence);
	fm->setf0Bounds(ssfpSequence->bandwidth());
	apply->SetAlgorithm(fm);
	apply->SetDataInput(0, ssfpSlice->GetOutput());
	apply->SetConstInput(0, T1Slice->GetOutput());
	if (B1) {
		B1Slice->SetInput(B1->GetOutput());
		B1Slice->SetRegionOfInterest(region);
		apply->SetConstInput(1, B1Slice->GetOutput());
	}
	if (mask) {
		maskSlice->SetInput(mask->GetOutput());
		maskSlice->SetRegionOfInterest(region);
		apply->SetMask(maskSlice->GetOutput());
	}
	time_t startTime;
	if (verbose) startTime = QI::printStartTime();
	apply->Update();
	QI::printElapsedTime(startTime);
	if (verbose) cout << "Writing output files." << endl;
	outPrefix = outPrefix + "FM_";
	QI::writeResult(apply->GetOutput(0), outPrefix + "PD.nii");
	QI::writeResult(apply->GetOutput(1), outPrefix + "T2.nii");
	QI::writeResult(apply->GetOutput(2), outPrefix + "f0.nii");
	QI::writeResiduals(apply->GetResidOutput(), outPrefix, all_residuals);
	return EXIT_SUCCESS;
}
