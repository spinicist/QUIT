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
#include "Filters/ApplyAlgorithmSliceBySliceFilter.h"
#include "Model.h"
#include "Sequence.h"
#include "RegionContraction.h"


using namespace std;
using namespace Eigen;

class FMFunctor : public DenseFunctor<double> {
public:
	const shared_ptr<SequenceBase> m_sequence;
	shared_ptr<SCD> m_model;
	ArrayXd m_data;
	const double m_T1, m_B1;

	FMFunctor(const shared_ptr<SCD> m, const shared_ptr<SequenceBase> s, const ArrayXd &d, const double T1, const double B1) :
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

class FixT2 : public DenseFunctor<double> {
public:
	const shared_ptr<SequenceBase> m_sequence;
	shared_ptr<SCD> m_model;
	ArrayXd m_data;
	const double m_T1, m_B1;
	double m_T2;

	FixT2(const shared_ptr<SCD> m, const shared_ptr<SequenceBase> s, const ArrayXd &d, const double T1, const double T2, const double B1) :
		DenseFunctor<double>(2, s->size()),
		m_model(m), m_sequence(s), m_data(d),
		m_T1(T1), m_T2(T2), m_B1(B1)
	{
		assert(static_cast<size_t>(m_data.rows()) == values());
	}

	void setT2(double T2) { m_T2 = T2; }
	int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
		eigen_assert(diffs.size() == values());

		ArrayXd fullparams(5);
		fullparams << params(0), m_T1, m_T2, params(1), m_B1;
		ArrayXcd s = m_sequence->signal(m_model, fullparams);
		diffs = s.abs() - m_data;
		return 0;
	}
};

template<typename T>
class FMAlgo : public Algorithm<T> {
protected:
	const shared_ptr<SCD> m_model = make_shared<SCD>();
	shared_ptr<SSFPSimple> m_sequence;

public:
	typedef typename Algorithm<T>::TArray TArray;

	void setSequence(shared_ptr<SSFPSimple> s) { m_sequence = s; }

	size_t numInputs() const override  { return m_sequence->count(); }
	size_t numConsts() const override  { return 2; }
	size_t numOutputs() const override { return 3; }
	size_t dataSize() const override   { return m_sequence->size(); }

	virtual TArray defaultConsts() {
		// T1 & B1
		TArray def = TArray::Ones(2);
		return def;
	}
};

class FMLMAlgo : public FMAlgo<double> {
protected:
	static const int m_iterations = 15;

public:
	virtual void apply(const TInput &data, const TArray &inputs,
					   TArray &outputs, TArray &resids) const override
	{
		const double T1 = inputs[0];
		if (isfinite(T1) && (T1 > 0.001)) {
			const double B1 = inputs[1];

            double bestF = numeric_limits<double>::infinity();
			for (int j = 1; j < 3; j++) {
				const double T2 = 0.045 * T1 * j; // From a Yarnykh paper T2/T1 = 0.045 in brain at 3T. Try the longer value for CSF
                Array2d b = m_sequence->bandwidth();
                double bw = b[1] - b[0];
                for (float f0guess = b[0] / 10.; f0guess < b[1] / 10.; f0guess += bw / 20.) {
					// First fix T2 and f0 to different starting points
					FixT2 fixT2(m_model, m_sequence, data, T1, T2, B1);
					NumericalDiff<FixT2> fixT2Diff(fixT2);
					LevenbergMarquardt<NumericalDiff<FixT2>> fixT2LM(fixT2Diff);
					fixT2LM.setMaxfev(m_iterations * (m_sequence->size() + 1));
                    VectorXd g(2); g << data.maxCoeff() * 2.5, f0guess;
					fixT2LM.minimize(g);

					// Now fit everything together
					FMFunctor full(m_model, m_sequence, data, T1, B1);
					NumericalDiff<FMFunctor> fullDiff(full);
					LevenbergMarquardt<NumericalDiff<FMFunctor>> fullLM(fullDiff);
					VectorXd fullP(3); fullP << g[0], T2, g[1]; // Now include T2
					fullLM.minimize(fullP);

                    double F = fullLM.fnorm();
					if (F < bestF) {
						outputs = fullP;
						bestF = F;
					}
				}
			}
			// PD sometimes go -ve, which is perfectly valid from the maths
			if (outputs[0] < 0)
				outputs[0] = -outputs[0];
			outputs[1] = clamp(outputs[1], 0.001, T1);
			VectorXd pfull(5); pfull << outputs[0], T1, outputs[1], outputs[2], B1; // Now include EVERYTHING to get a residual
			ArrayXd theory = m_sequence->signal(m_model, pfull).abs();
			resids = (data - theory);
		} else {
			// No point in processing -ve T1
			outputs.setZero();
			resids.setZero();
		}
	}
};

class XFMFunctor : public DenseFunctor<double> {
public:
	const shared_ptr<SequenceBase> m_sequence;
	shared_ptr<SCD> m_model;
	ArrayXcd m_data;
	const double m_T1, m_B1;

	XFMFunctor(const shared_ptr<SCD> m, const shared_ptr<SequenceBase> s, const ArrayXcd &d, const double T1, const double B1) :
		DenseFunctor<double>(3, s->size()),
		m_model(m), m_sequence(s), m_data(d),
		m_T1(T1), m_B1(B1)
	{
		assert(static_cast<size_t>(m_data.rows()) == values());
	}

	int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
		eigen_assert(diffs.size() == values());

		ArrayXd fullparams(5);
		fullparams << params(0), m_T1, params(1), params(2), m_B1;
		ArrayXcd s = m_sequence->signal(m_model, fullparams);
		diffs = (s - m_data).abs();
		return 0;
	}
};

class XFMLMAlgo : public FMAlgo<complex<double>> {
protected:
	static const int m_iterations = 15;

public:
	virtual void apply(const TInput &data, const TArray &inputs,
					   TArray &op, TArray &resids) const override
	{
		const double T1 = inputs[0];
		const double B1 = inputs[1];
		if (isfinite(T1) && (T1 > 0.001)) {
			complex<double> avg = data.mean();
			XFMFunctor full(m_model, m_sequence, data, T1, B1);
			NumericalDiff<XFMFunctor> fullDiff(full);
            LevenbergMarquardt<NumericalDiff<XFMFunctor>> fullLM(fullDiff);
            double bestF = numeric_limits<double>::infinity();
            VectorXd bestP;
            double bw = (2. / m_sequence->TR());
            for (int f = -1; f < 2; f++) {
                double f0guess = arg(avg) / (2. * M_PI * m_sequence->TR()) + f * bw;
                VectorXd P3(3); P3 << abs(avg) * 5.0, 0.045 * T1, f0guess;
                fullLM.minimize(P3);
                if (fullLM.fnorm() < bestF) {
                    bestP = P3;
                    bestF = fullLM.fnorm();
                }
            }
            //P[1] = clamp(P[1], 0.001, T1);
            VectorXd pfull(5); pfull << bestP[0], T1, bestP[1], bestP[2], B1; // Now include EVERYTHING to get a residual
			ArrayXcd theory = m_sequence->signal(m_model, pfull);
			resids = (data - theory).abs();
            op = bestP;
		} else {
			// No point in processing -ve T1
			op.setZero();
			resids.setZero();
		}
	}
};

class FMSRCAlgo : public FMAlgo<double> {
private:
	size_t m_samples = 2000, m_retain = 20, m_contractions = 10;
	Array2d m_f0Bounds = Array2d::Zero();

public:
	void setRCPars(size_t c, size_t s, size_t r) { m_contractions = c; m_samples = s; m_retain = r; }
	void setf0Bounds(Array2d b) { m_f0Bounds = b; }

	virtual void apply(const TInput &data, const TArray &inputs,
					   TArray &outputs, TArray &resids) const override
	{
		ArrayXd thresh(3); thresh.setConstant(0.05);
		ArrayXd weights(m_sequence->size()); weights.setOnes();
		ArrayXXd bounds = ArrayXXd::Zero(3, 2);
		double T1 = inputs[0];
		if (isfinite(T1) && (T1 > 0.001)) {
			double B1 = inputs[1];
			bounds(0, 0) = 0.;
			bounds(0, 1) = data.array().abs().maxCoeff() * 25;
			bounds(1,0) = 0.001;
			bounds(1,1) = T1;
			bounds.row(2) = m_sequence->bandwidth();
			//cout << "T1 " << T1 << " B1 " << B1 << " inputs " << inputs.transpose() << endl;
			//cout << bounds << endl;
			FMFunctor func(m_model, m_sequence, data, T1, B1);
			RegionContraction<FMFunctor> rc(func, bounds, weights, thresh,
											m_samples, m_retain, m_contractions, 0.02, true, false);
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

	int start_slice = 0, stop_slice = 0;
	int verbose = false, prompt = true, all_residuals = false,
	    fitFinite = false, flipData = false, use_complex = false, use_src = false;
	string outPrefix;
	QI::ReadImageF::Pointer mask = ITK_NULLPTR, B1 = ITK_NULLPTR;

	const string usage {
	"Usage is: despot2-fm [options] T1_map ssfp_file\n\
	\
	Options:\n\
		--help, -h        : Print this message\n\
		--verbose, -v     : Print slice processing times\n\
		--no-prompt, -n   : Suppress input prompts\n\
		--mask, -m file   : Mask input with specified file\n\
		--out, -o path    : Add a prefix to the output filenames\n\
		--B1, -b file     : B1 Map file (ratio)\n\
		--algo, -a l      : Use 2-step LM algorithm (default)\n\
				   s      : Use Stochastic Region Contraction\n\
				   x      : Fit to complex data\n\
		--start, -s N     : Start processing from slice N\n\
		--stop, -p  N     : Stop processing at slice N\n\
		--flip, -F        : Data order is phase, then flip-angle (default opposite)\n\
		--finite          : Use finite pulse length correction\n\
		--resids, -r      : Write out per flip-angle residuals\n\
		--threads, -T N   : Use N threads (default=hardware limit)\n"
	};
	/* --complex, -x     : Fit to complex data\n\ */

	struct option long_opts[] = {
		{"help", no_argument, 0, 'h'},
		{"verbose", no_argument, 0, 'v'},
		{"no-prompt", no_argument, 0, 'n'},
		{"mask", required_argument, 0, 'm'},
		{"out", required_argument, 0, 'o'},
		{"B1", required_argument, 0, 'b'},
		{"algo", required_argument, 0, 'a'},
		{"start", required_argument, 0, 's'},
		{"stop", required_argument, 0, 'p'},
		{"flip", required_argument, 0, 'F'},
		{"threads", required_argument, 0, 'T'},
		{"finite", no_argument, &fitFinite, 1},
		{"resids", no_argument, 0, 'r'},
		{0, 0, 0, 0}
	};
	const char* short_opts = "hvnm:o:b:a:s:p:FT:rd:";

	// Do first pass to get the algorithm type, then do everything else
	int indexptr = 0;
	char c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'a':
			switch (*optarg) {
				case 'l': if (verbose) cout << "LM algorithm selected." << endl; break;
				case 's': use_src = true; if (verbose) cout << "SRC algorithm selected." << endl; break;
				case 'x':
					use_complex = true;
					if (verbose) cout << "XFM algorithm selected." << endl;
					break;
				default: throw(runtime_error(string("Unknown algorithm type ") + optarg)); break;
			} break;
			default: break;
		}
	}

	optind = 1;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
			case 'v': case 'n': case 'a': break; // Already handled
			case 'm':
				if (verbose) cout << "Reading mask file " << optarg << endl;
				mask = QI::ReadImageF::New();
				mask->SetFileName(optarg);
				break;
			case 'o':
				outPrefix = optarg;
				if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'b':
				if (verbose) cout << "Reading B1 file: " << optarg << endl;
				B1 = QI::ReadImageF::New();
				B1->SetFileName(optarg);
				break;
			case 's': start_slice = atoi(optarg); break;
			case 'p': stop_slice = atoi(optarg); break;
			case 'F': flipData = true; break;
			case 'T': itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg)); break;
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
	if ((argc - optind) != 2) {
		cout << usage << endl;
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
	auto T1 = QI::ReadImageF::New();
	T1->SetFileName(argv[optind++]);

	if (use_complex) {
        if (verbose) {
            cout << "Running complex pipelin." << endl;
            cout << "Opening SSFP file: " << argv[optind] << endl;
        }
		auto ssfpFile = QI::ReadTimeseriesXF::New();
		auto ssfpData = QI::TimeseriesToVectorXF::New();
		auto ssfpFlip = QI::ReorderXF::New();
		ssfpFile->SetFileName(argv[optind++]);
		ssfpData->SetInput(ssfpFile->GetOutput());
		ssfpFlip->SetInput(ssfpData->GetOutput());
		if (flipData) {
			ssfpFlip->SetStride(ssfpSequence->phases());
		}
		typedef itk::ApplyAlgorithmSliceBySliceFilter<XFMLMAlgo, complex<float>, float, 3> TFMFilter;
		auto apply = TFMFilter::New();
		shared_ptr<XFMLMAlgo> fm = make_shared<XFMLMAlgo>();
		fm->setSequence(ssfpSequence);
		apply->SetAlgorithm(fm);
		apply->SetInput(0, ssfpFlip->GetOutput());
		apply->SetConst(0, T1->GetOutput());
		apply->SetSlices(start_slice, stop_slice);
		if (B1) {
			apply->SetConst(1, B1->GetOutput());
		}
		if (mask) {
			apply->SetMask(mask->GetOutput());
		}
		time_t startTime;
		if (verbose) {
			startTime = QI::printStartTime();
			auto monitor = QI::SliceMonitor<TFMFilter>::New();
			apply->AddObserver(itk::IterationEvent(), monitor);
		}
		apply->Update();
		if (verbose) {
			QI::printElapsedTime(startTime);
            cout << "Writing output files. Prefix is " << outPrefix << endl;
		}
		outPrefix = outPrefix + "FM_";
		QI::writeResult(apply->GetOutput(0), outPrefix + "PD.nii");
		QI::writeResult(apply->GetOutput(1), outPrefix + "T2.nii");
		QI::writeResult(apply->GetOutput(2), outPrefix + "f0.nii");
		QI::writeResiduals(apply->GetResidOutput(), outPrefix, all_residuals);
	} else {
		if (verbose) cout << "Opening SSFP file: " << argv[optind] << endl;
		auto ssfpFile = QI::ReadTimeseriesF::New();
		auto ssfpData = QI::TimeseriesToVectorF::New();
		auto ssfpFlip = QI::ReorderF::New();
		ssfpFile->SetFileName(argv[optind++]);
		ssfpData->SetInput(ssfpFile->GetOutput());
		ssfpFlip->SetInput(ssfpData->GetOutput());
		if (flipData) {
			ssfpFlip->SetStride(ssfpSequence->phases());
		}
		typedef itk::ApplyAlgorithmSliceBySliceFilter<FMAlgo<double>, float, float, 3> TFMFilter;
		auto apply = TFMFilter::New();
		shared_ptr<FMAlgo<double>> fm;
		if (use_src)
			fm = make_shared<FMSRCAlgo>();
		else
			fm = make_shared<FMLMAlgo>();
		fm->setSequence(ssfpSequence);
		apply->SetAlgorithm(fm);
		apply->SetInput(0, ssfpFlip->GetOutput());
		apply->SetConst(0, T1->GetOutput());
		apply->SetSlices(start_slice, stop_slice);
		if (B1) {
			apply->SetConst(1, B1->GetOutput());
		}
		if (mask) {
			apply->SetMask(mask->GetOutput());
		}
		time_t startTime;
		if (verbose) {
			startTime = QI::printStartTime();
			auto monitor = QI::SliceMonitor<TFMFilter>::New();
			apply->AddObserver(itk::IterationEvent(), monitor);
		}
		apply->Update();
		if (verbose) {
			QI::printElapsedTime(startTime);
            cout << "Writing output files. Prefix is " << outPrefix << endl;
		}
		outPrefix = outPrefix + "FM_";
		QI::writeResult(apply->GetOutput(0), outPrefix + "PD.nii");
		QI::writeResult(apply->GetOutput(1), outPrefix + "T2.nii");
		QI::writeResult(apply->GetOutput(2), outPrefix + "f0.nii");
		QI::writeResiduals(apply->GetResidOutput(), outPrefix, all_residuals);
	}
	return EXIT_SUCCESS;
}
