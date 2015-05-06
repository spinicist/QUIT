/*
 *  despot1.cpp
 *  Part of quitk
 *
 *  Created by Tobias Wood on 12/01/2015.
 *  Copyright (c) 2011-2013 Tobias Wood.
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
#include "itkImageFileWriter.h"
#include "itkVectorImage.h"
#include "itkImageToImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkComposeImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "Filters/ImageToVectorFilter.h"
#include "Filters/ApplyAlgorithmFilter.h"
#include "Model.h"
#include "Sequence.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
// Algorithm Subclass
//******************************************************************************
class DESPOT1 : public Algorithm<double> {
	public:
		enum class Type { LLS, WLLS, NLLS };
	private:
		Type m_type = Type::LLS;
		size_t m_iterations = 4;

		// T1 only Functor
		class T1Functor : public DenseFunctor<double> {
			protected:
				const shared_ptr<SequenceBase> m_sequence;
				const ArrayXd m_data;
				const bool m_debug;
				const double m_B1;
				const shared_ptr<SCD> m_model = make_shared<SCD>();

			public:
				T1Functor(const shared_ptr<SequenceBase> cs, const ArrayXd &data,
						  const double B1, const bool debug) :
					DenseFunctor<double>(2, cs->size()),
					m_sequence(cs), m_data(data),
					m_B1(B1), m_debug(debug)
				{
					assert(static_cast<size_t>(m_data.rows()) == values());
				}

				int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
					eigen_assert(diffs.size() == values());
					VectorXd fullParams = VectorXd::Zero(5);
					fullParams.head(2) = params;
					fullParams(4) = m_B1;
					ArrayXcd s = m_sequence->signal(m_model, fullParams);
					diffs = s.abs() - m_data;
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
		void setType(Type t) { m_type = t; }
		void setIterations(size_t n) { m_iterations = n; }

		size_t numConsts() const override { return 1; }
		size_t numOutputs() const override { return 2; }

		virtual VectorXd defaultConsts() {
			// Only B1, default is 1
			VectorXd def = VectorXd::Ones(1);
			return def;
		}

		virtual void apply(const shared_ptr<SequenceBase> sequence,
		                   const VectorXd &data,
		                   const VectorXd &inputs,
		                   VectorXd &outputs,
		                   ArrayXd &resids) const override
		{
			double B1 = inputs[0];
			ArrayXd flip = sequence->flip() * B1;
			VectorXd Y = data.array() / flip.sin();
			MatrixXd X(Y.rows(), 2);
			X.col(0) = data.array() / flip.tan();
			X.col(1).setOnes();
			VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
			outputs[1] = -sequence->TR() / log(b[0]);
			outputs[0] = b[1] / (1. - b[0]);
			if (m_type == Type::WLLS) {
				VectorXd W(sequence->size());
				for (size_t n = 0; n < m_iterations; n++) {
					W = (flip.sin() / (1. - (exp(-sequence->TR()/outputs[1])*flip.cos()))).square();
					b = (X.transpose() * W.asDiagonal() * X).partialPivLu().solve(X.transpose() * W.asDiagonal() * Y);
					outputs[1] = -sequence->TR() / log(b[0]);
					outputs[0] = b[1] / (1. - b[0]);
				}
			} else if (m_type == Type::NLLS) {
				T1Functor f(sequence, data, B1, false);
				NumericalDiff<T1Functor> nDiff(f);
				LevenbergMarquardt<NumericalDiff<T1Functor>> lm(nDiff);
				lm.setMaxfev(m_iterations * (sequence->size() + 1));
				lm.minimize(outputs);
			}

			ArrayXd theory = One_SPGR(sequence->flip(), sequence->TR(), outputs[0], outputs[1], B1).array().abs();
			resids = (data.array() - theory);
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
	--its, -i N       : Max iterations for WLLS (default 4)\n\
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
	{"threads", required_argument, 0, 'T'},
	{"resids", no_argument, 0, 'r'},
	{0, 0, 0, 0}
};
static const char *short_opts = "hvnm:o:b:a:i:T:r";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	try { // To fix uncaught exceptions on Mac
	//cout << version << endl << credit_shared << endl;
	Eigen::initParallel();

	typedef itk::Image<float, 3> FloatImage;
	typedef itk::VectorImage<float, 3> FloatVectorImage;

	typedef itk::ImageFileReader<FloatImage> Reader;
	typedef itk::ImageFileReader<itk::Image<float, 4>> Reader4D;
	typedef itk::ImageFileWriter<FloatImage> Writer;

	Reader::Pointer mask = ITK_NULLPTR;
	Reader::Pointer B1   = ITK_NULLPTR;

	shared_ptr<DESPOT1> algo = make_shared<DESPOT1>();

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'm':
				cout << "Opening mask file " << optarg << endl;
				mask = Reader::New();
				mask->SetFileName(optarg);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'b':
				cout << "Opening B1 file: " << optarg << endl;
				B1 = Reader::New();
				B1->SetFileName(optarg);
				break;
			case 'a':
				switch (*optarg) {
					case 'l': algo->setType(DESPOT1::Type::LLS);  cout << "LLS algorithm selected." << endl; break;
					case 'w': algo->setType(DESPOT1::Type::WLLS); cout << "WLLS algorithm selected." << endl; break;
					case 'n': algo->setType(DESPOT1::Type::NLLS); cout << "NLLS algorithm selected." << endl; break;
					default:
						cout << "Unknown algorithm type " << optarg << endl;
						return EXIT_FAILURE;
						break;
				} break;
			case 'i':
				nIterations = atoi(optarg);
				break;
			case 'T':
				itk::MultiThreader::SetGlobalDefaultNumberOfThreads(atoi(optarg));
				break;
			case 'r': all_residuals = true; break;
			case 'h':
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 1) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		return EXIT_FAILURE;
	}

	if (verbose) cout << "Opening SPGR file: " << argv[optind] << endl;
	// Deal with possible complex data
	string inputFilename = argv[optind++];
	auto input = itk::ImageFileReader<itk::Image<float, 4>>::New();
	auto cinput = itk::ImageFileReader<itk::Image<complex<float>, 4>>::New();
	auto mag = itk::ComplexToModulusImageFilter<itk::Image<complex<float>, 4>, itk::Image<float, 4>>::New();
	auto convert = ImageToVectorFilter<float>::New();
	auto imageIO = itk::ImageIOFactory::CreateImageIO(inputFilename.c_str(), itk::ImageIOFactory::ReadMode);
	imageIO->SetFileName(inputFilename);
	imageIO->ReadImageInformation();
	const itk::ImageIOBase::IOPixelType pType = imageIO->GetPixelType();
	switch (pType) {
		case itk::ImageIOBase::SCALAR: {
			// Just read as is
			input->SetFileName(inputFilename);
			convert->SetInput(input->GetOutput());
		} break;
		case itk::ImageIOBase::COMPLEX: {
			// Convert to magnitude first
			cinput->SetFileName(inputFilename);
			mag->SetInput(cinput->GetOutput());
			convert->SetInput(mag->GetOutput());
		} break;
		default:
			throw(runtime_error("Pixel Type (" + imageIO->GetPixelTypeAsString(pType) + ") not supported."));
			break;
	}

	shared_ptr<SPGRSimple> spgrSequence = make_shared<SPGRSimple>(prompt);
	if (verbose) cout << *spgrSequence;

	auto d1 = itk::ApplyAlgorithmFilter<float, DESPOT1>::New();
	d1->SetSequence(spgrSequence);
	d1->SetAlgorithm(algo);
	d1->Setup();
	d1->SetDataInput(0, convert->GetOutput());
	if (mask) {
		cout << "Setting mask" << endl;
		d1->SetMask(mask->GetOutput());
	}
	if (B1)
		d1->SetConstInput(0, B1->GetOutput());
	algo->setIterations(nIterations);
	if (verbose) cout << "Processing" << endl;
	d1->Update();

	if (verbose)
		cout << "Writing results." << endl;
	outPrefix = outPrefix + "D1_";

	Writer::Pointer T1File = Writer::New();
	Writer::Pointer PDFile = Writer::New();
	Writer::Pointer ResFile = Writer::New();

	T1File->SetFileName(outPrefix + "T1.nii");
	PDFile->SetFileName(outPrefix + "PD.nii");
	ResFile->SetFileName(outPrefix + "residual.nii");

	PDFile->SetInput(d1->GetOutput(0));
	T1File->SetInput(d1->GetOutput(1));

	auto magFilter = itk::VectorMagnitudeImageFilter<FloatVectorImage, FloatImage>::New();
	magFilter->SetInput(d1->GetResidOutput());
	ResFile->SetInput(magFilter->GetOutput());

	T1File->Update();
	PDFile->Update();
	ResFile->Update();

	if (all_residuals) {
		auto vecWriter = itk::ImageFileWriter<FloatVectorImage>::New();
		vecWriter->SetFileName(outPrefix + "residuals.nii");
		vecWriter->SetInput(d1->GetResidOutput());
		vecWriter->Update();
	}
	cout << "All done." << endl;
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
