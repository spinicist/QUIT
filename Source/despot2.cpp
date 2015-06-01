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

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImage.h"
#include "itkImageToImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkComposeImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "Filters/ImageToVectorFilter.h"
#include "Filters/ApplyAlgorithmFilter.h"
#include "Model.h"
#include "Sequence.h"
#include "Util.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
// Algorithm Subclass
//******************************************************************************
class D2 : public Algorithm<complex<double>> {
	public:
		enum class Type { LLS, WLLS, NLLS };
	private:
		Type m_type = Type::LLS;
		size_t m_iterations = 10;
		bool m_elliptical = false, m_negflip;
		double m_thresh = -numeric_limits<double>::infinity();
		double m_lo = -numeric_limits<double>::infinity();
		double m_hi = numeric_limits<double>::infinity();

		//******************************************************************************
		// T2 Only Functor
		//******************************************************************************
		class D2Functor : public DenseFunctor<double> {
			public:
				const shared_ptr<SequenceBase> m_sequence;
				ArrayXcd m_data;
				const double m_T1, m_B1;
				const bool m_complex, m_debug;
				const shared_ptr<SCD> m_model = make_shared<SCD>();

				D2Functor(const double T1, const shared_ptr<SequenceBase> s, const ArrayXcd &d, const double B1, const bool fitComplex, const bool debug = false) :
					DenseFunctor<double>(3, s->size()),
					m_sequence(s), m_data(d), m_complex(fitComplex), m_debug(debug),
					m_T1(T1), m_B1(B1)
				{
					assert(static_cast<size_t>(m_data.rows()) == values());
				}

				int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
					eigen_assert(diffs.size() == values());

					ArrayXd fullparams(5);
					fullparams << params(0), m_T1, params(1), params(2), m_B1;
					ArrayXcd s = m_sequence->signal(m_model, fullparams);
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

	public:
		void setElliptical(bool e) { m_elliptical = e; }
		void setType(Type t) { m_type = t; }
		void setIterations(size_t n) { m_iterations = n; }
		void setThreshold(double t) { m_thresh = t; }
		void setClamp(double lo, double hi) { m_lo = lo; m_hi = hi; }

		size_t numConsts() const override { return 2; } // T1, B1
		size_t numOutputs() const override { return 3; } // PD, T2, offRes (for elliptical)

		virtual VectorXd defaultConsts() {
			//T1 & B1
			VectorXd def = VectorXd::Ones(2);
			return def;
		}

		virtual void apply(const shared_ptr<SequenceBase> sequence,
		                   const VectorXcd &data,
		                   const VectorXd &constants,
		                   VectorXd &outputs,
		                   ArrayXd &resids) const override
		{
			double T2, E2, PD, offRes;
			const double TR = sequence->TR();
			const double T1 = constants[0];
			const double B1 = constants[1];
			const double E1 = exp(-TR / T1);
			const ArrayXd localAngles = (sequence->flip() * B1);
			const ArrayXd s = data.array().abs();
			VectorXd Y = s / localAngles.sin();
			MatrixXd X(Y.rows(), 2);
			X.col(0) = s / localAngles.tan();
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
			if (m_type == Type::WLLS) {
				VectorXd W(sequence->size());
				for (size_t n = 0; n < m_iterations; n++) {
					if (m_elliptical) {
						W = ((1. - E1*E2) * localAngles.sin() / (1. - E1*E2*E2 - (E1 - E2*E2)*localAngles.cos())).square();
					} else {
						W = ((1. - E1*E2) * localAngles.sin() / (1. - E1*E2 - (E1 - E2)*localAngles.cos())).square();
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
			} else if (m_type == Type::NLLS) {
				D2Functor f(T1, sequence, data, B1, false, false);
				NumericalDiff<D2Functor> nDiff(f);
				LevenbergMarquardt<NumericalDiff<D2Functor>> lm(nDiff);
				lm.setMaxfev(m_iterations * (sequence->size() + 1));
				VectorXd p(3);
				p << PD, T2, offRes;
				lm.minimize(p);
				PD = p(0); T2 = p(1); offRes = p(2);
			}
			if (PD < m_thresh) {
				PD = 0.;
				T2 = 0.;
				offRes = 0.;
			}
			T2 = clamp(T2, m_lo, m_hi);

			if (m_elliptical) {
				// Use phase of mean instead of mean of phase to avoid wrap issues
				if (m_negflip)
					offRes = -arg(data.mean()) / (M_PI * TR);
				else
					offRes = arg(data.mean() / (M_PI * TR));
			}

			outputs[0] = PD;
			outputs[1] = T2;
			outputs[3] = offRes;

			VectorXd p(5); p << PD, T1, T2, offRes, B1;
			shared_ptr<SCD> model = make_shared<SCD>();
			ArrayXd theory = sequence->signal(model, p).abs();
			resids = (s - theory);

		}
};


//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: despot2 [options] T1_map ssfp_file\n\
\n\
Options:\n\
	--help, -h        : Print this message.\n\
	--verbose, -v     : Print slice processing times.\n\
	--no-prompt, -n   : Suppress input prompts.\n\
	--mask, -m file   : Mask input with specified file.\n\
	--out, -o path    : Add a prefix to the output filenames.\n\
	--B1 file         : B1 Map file.\n\
	--thresh, -t n    : Threshold maps at PD < n\n\
	--clamp, -c n     : Clamp T1 between 0 and n\n\
	--algo, -a l      : LLS algorithm (default)\n\
	           w      : WLLS algorithm\n\
	           n      : NLLS (Levenberg-Marquardt)\n\
	--its, -i N       : Max iterations for WLLS / NLLS (default 10)\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n\
	--elliptical, -e  : Input is band-free elliptical data.\n\
	Only applies for elliptical data:\n\
	--negflip         : Flip-angles are in negative sense.\n"
};

enum class Algos { LLS, WLLS, NLLS };
static int verbose = false, prompt = true, elliptical = false, negflip = false, all_residuals = false;
static Algos algo;
static string outPrefix;
static struct option long_opts[] =
{
	{"B1", required_argument, 0, '1'},
	{"elliptical", no_argument, 0, 'e'},
	{"negflip", no_argument, &negflip, true},
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
static const char *short_opts = "hm:o:b:t:c:vna:i:T:er";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	try { // To fix uncaught exceptions on Mac

	Eigen::initParallel();

	typedef itk::Image<float, 3> FloatImage;
	typedef itk::VectorImage<float, 3> FloatVectorImage;
	typedef itk::Image<complex<float>, 3> CFloatImage;
	typedef itk::VectorImage<complex<float>, 3> CFloatVectorImage;
	typedef itk::ImageFileReader<FloatImage> Reader;

	typedef itk::ImageFileWriter<FloatImage> Writer;

	Reader::Pointer mask = ITK_NULLPTR;
	Reader::Pointer B1   = ITK_NULLPTR;

	shared_ptr<D2> algo = make_shared<D2>();

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
			case 'o':
				outPrefix = optarg;
				if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'm':
				if (verbose) cout << "Reading mask file " << optarg << endl;
				mask = Reader::New();
				mask->SetFileName(optarg);
				break;
			case 'b':
				if (verbose) cout << "Reading B1 file: " << optarg << endl;
				B1 = Reader::New();
				B1->SetFileName(optarg);
				break;
			case 't': algo->setThreshold(atof(optarg)); break;
			case 'c': algo->setClamp(0, atof(optarg)); break;
			case 'a':
				switch (*optarg) {
					case 'l': algo->setType(D2::Type::LLS);  cout << "LLS algorithm selected." << endl; break;
					case 'w': algo->setType(D2::Type::WLLS); cout << "WLLS algorithm selected." << endl; break;
					case 'n': algo->setType(D2::Type::NLLS); cout << "NLLS algorithm selected." << endl; break;
					default:
						cout << "Unknown algorithm type " << optarg << endl;
						return EXIT_FAILURE;
						break;
				} break;
			case 'i':
				algo->setIterations(atoi(optarg));
				break;
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'e': elliptical = true; algo->setElliptical(elliptical); break;
			case 'r': all_residuals = true; break;
			case 'T':
				itk::MultiThreader::SetGlobalDefaultNumberOfThreads(atoi(optarg));
				break;
			case 0:
				// Just a flag
				break;
			case '?': // getopt will print an error message
			case 'h':
				cout << usage << endl;
				return EXIT_SUCCESS;
		}
	}
	//if (verbose) cout << version << endl << credit_shared << endl;
	if ((argc - optind) != 2) {
		cout << "Wrong number of arguments. Need a T1 map and SSFP file." << endl;
		cout << usage << endl;
		return EXIT_FAILURE;
	}
	if (verbose) cout << "Reading T1 Map from: " << argv[optind] << endl;
	Reader::Pointer T1File = Reader::New();
	T1File->SetFileName(argv[optind++]);

	if (verbose) cout << "Opening SSFP file: " << argv[optind] << endl;
	auto ssfp4D = QUITK::ReadXFloatTimeseries::New();
	ssfp4D->SetFileName(argv[optind++]);
	auto ssfp3D = itk::ImageToVectorFilter<QUITK::XFloatTimeseries>::New();
	ssfp3D->SetInput(ssfp4D->GetOutput());

	shared_ptr<SteadyState> ssfp;
	if (elliptical) {
		ssfp = make_shared<SSFPEllipse>(prompt);
	} else {
		ssfp = make_shared<SSFPSimple>(prompt);
	}
	if (verbose) cout << *ssfp << endl;

	auto DESPOT2 = itk::ApplyAlgorithmFilter<complex<float>, D2>::New();
	DESPOT2->SetSequence(ssfp);
	DESPOT2->SetAlgorithm(algo);
	DESPOT2->Setup();
	DESPOT2->SetDataInput(0, ssfp3D->GetOutput());
	DESPOT2->SetConstInput(0, T1File->GetOutput());
	if (B1)
		DESPOT2->SetConstInput(1, B1->GetOutput());
	if (mask)
		DESPOT2->SetMask(mask->GetOutput());

	time_t startTime;
	if (verbose) {
		cout << "DESPOT2 setup complete. Processing." << endl;
		startTime = QUITK::printStartTime();
	}
	DESPOT2->Update();
	if (verbose) {
		QUITK::printElapsedTime(startTime);
		cout << "Writing results." << endl;
	}

	outPrefix = outPrefix + "D2_";

	Writer::Pointer PDFile = Writer::New();
	Writer::Pointer T2File = Writer::New();
	Writer::Pointer f0File = Writer::New();
	Writer::Pointer ResFile = Writer::New();

	PDFile->SetFileName(outPrefix + "PD.nii");
	PDFile->SetInput(DESPOT2->GetOutput(0));
	PDFile->Update();
	T2File->SetFileName(outPrefix + "T2.nii");
	T2File->SetInput(DESPOT2->GetOutput(1));
	T2File->Update();
	if (elliptical) {
		f0File->SetFileName(outPrefix + "f0.nii");
		f0File->SetInput(DESPOT2->GetOutput(2));
		f0File->Update();
	}

	auto magFilter = itk::VectorMagnitudeImageFilter<FloatVectorImage, FloatImage>::New();
	magFilter->SetInput(DESPOT2->GetResidOutput());
	ResFile->SetFileName(outPrefix + "residual.nii");
	ResFile->SetInput(magFilter->GetOutput());
	ResFile->Update();

	if (all_residuals) {
		auto vecWriter = itk::ImageFileWriter<FloatVectorImage>::New();
		vecWriter->SetFileName(outPrefix + "residuals.nii");
		vecWriter->SetInput(DESPOT2->GetResidOutput());
		vecWriter->Update();
	}
	if (verbose) cout << "All done." << endl;
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
