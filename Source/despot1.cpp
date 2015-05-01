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

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImage.h"
#include "itkImageToImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkComposeImageFilter.h"

#include "Filters/ImageToVectorFilter.h"
#include "Model.h"
#include "Sequence.h"

using namespace std;
using namespace Eigen;

namespace itk{

template<typename TVectorImage, typename TImage>
class DESPOT1Filter : public ImageToImageFilter<TVectorImage, TImage>
{
public:
	/** Standard class typedefs. */
	typedef DESPOT1Filter                      Self;
	typedef ImageToImageFilter<TVectorImage, TImage> Superclass;
	typedef SmartPointer<Self>                 Pointer;

	enum class Algos { LLS, WLLS, NLLS };

	itkNewMacro(Self); /** Method for creation through the object factory. */
	itkTypeMacro(DESPOT1Filter, ImageToImageFilter); /** Run-time type information (and related methods). */

	void SetInput(const TVectorImage *SPGR);
	typename TVectorImage::ConstPointer GetInput();
	void SetMask(const TImage *mask);
	typename TImage::ConstPointer GetMask();
	void SetB1(const TImage *B1);
	typename TImage::ConstPointer GetB1();
	TImage *GetOutputT1();
	TImage *GetOutputPD();
	TImage *GetOutputRes();

	void SetSequence(const SPGRSimple &seq);
	void SetAlgorithm(const Algos &a);
	void SetIterations(const size_t &n);

	virtual void Update();

	void PrintDirections() {
		cout << __PRETTY_FUNCTION__ << endl;
		typedef ImageBase< 3 > ImageBaseType;
		ImageBaseType *ptr = ITK_NULLPTR;
		InputDataObjectIterator it(this);

		for(; !it.IsAtEnd(); ++it ) {
			ptr = dynamic_cast< ImageBaseType * >( it.GetInput() );
			cout << it.GetName() << endl << ptr->GetDirection() << endl;
		}
	}

protected:
	DESPOT1Filter();
	~DESPOT1Filter(){}

	virtual void GenerateData(); // Does the work
	DataObject::Pointer MakeOutput(unsigned int idx); // Create the Output

private:
	DESPOT1Filter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented

	SPGRSimple m_sequence;
	Algos m_algorithm;
	size_t m_iterations;
};
}

#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk {



template<typename TVectorImage, typename TImage>
DESPOT1Filter<TVectorImage, TImage>::DESPOT1Filter() :
	m_sequence({},0)
{
	cout << __PRETTY_FUNCTION__ << endl;
	this->SetNumberOfRequiredInputs(3);
	this->SetNumberOfRequiredOutputs(3);

	this->SetNthOutput(0, this->MakeOutput(0));
	this->SetNthOutput(1, this->MakeOutput(1));
	this->SetNthOutput(2, this->MakeOutput(2));
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetInput(const TVectorImage *image) {
	cout << __PRETTY_FUNCTION__ << endl;
	this->SetNthInput(0, const_cast<TVectorImage*>(image));
	cout << this->GetInput()->GetDirection() << endl;
}
template< typename TVectorImage, typename TImage>
typename TVectorImage::ConstPointer DESPOT1Filter<TVectorImage, TImage>::GetInput() {
	cout << __PRETTY_FUNCTION__ << endl;
	return static_cast<const TVectorImage *> (this->ProcessObject::GetInput(0));
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetMask(const TImage* image) {
	cout << __PRETTY_FUNCTION__ << endl;
	this->SetNthInput(1, const_cast<TImage*>(image));
	cout << this->GetMask()->GetDirection() << endl;
}
template< typename TVectorImage, typename TImage>
typename TImage::ConstPointer DESPOT1Filter<TVectorImage, TImage>::GetMask() {
	cout << __PRETTY_FUNCTION__ << endl;
	return static_cast<const TImage *>(this->ProcessObject::GetInput(1));
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetB1(const TImage* image) {
	cout << __PRETTY_FUNCTION__ << endl;
	this->SetNthInput(2, const_cast<TImage*>(image));
	cout << this->GetB1()->GetDirection() << endl;
}
template<typename TVectorImage, typename TImage>
typename TImage::ConstPointer DESPOT1Filter<TVectorImage, TImage>::GetB1() {
	cout << __PRETTY_FUNCTION__ << endl;
	return static_cast<const TImage *>(this->ProcessObject::GetInput(2));
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetSequence(const SPGRSimple &seq) {
	cout << __PRETTY_FUNCTION__ << endl;
	m_sequence = seq;
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetAlgorithm(const Algos &a) {
	cout << __PRETTY_FUNCTION__ << endl;
	m_algorithm = a;
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetIterations(const size_t &n) {
	cout << __PRETTY_FUNCTION__ << endl;
	m_iterations = n;
}

template<typename TVectorImage, typename TImage>
DataObject::Pointer DESPOT1Filter<TVectorImage, TImage>::MakeOutput(unsigned int idx) {
	cout << __PRETTY_FUNCTION__ << endl;
	DataObject::Pointer output;

	switch ( idx ) {
	case 0: case 1: case 2: output = ( TImage::New() ).GetPointer(); break;
	default:
		std::cerr << "No output " << idx << std::endl;
		output = NULL;
		break;
	}
	return output.GetPointer();
}

template< typename TVectorImage, typename TImage>
TImage* DESPOT1Filter<TVectorImage, TImage>::GetOutputT1() {
	cout << __PRETTY_FUNCTION__ << endl;
	return dynamic_cast<TImage *>(this->ProcessObject::GetOutput(0) );
}
template< typename TVectorImage, typename TImage>
TImage* DESPOT1Filter<TVectorImage, TImage>::GetOutputPD() {
	cout << __PRETTY_FUNCTION__ << endl;
	return dynamic_cast<TImage *>(this->ProcessObject::GetOutput(1) );
}
template< typename TVectorImage, typename TImage>
TImage* DESPOT1Filter<TVectorImage, TImage>::GetOutputRes() {
	cout << __PRETTY_FUNCTION__ << endl;
	return dynamic_cast<TImage *>(this->ProcessObject::GetOutput(2) );
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::Update() {
	cout << __PRETTY_FUNCTION__ << endl;
	cout << static_cast<const TVectorImage *>(this->ProcessObject::GetInput(0))->GetDirection() << endl;
	cout << static_cast<const TImage *>(this->ProcessObject::GetInput(1))->GetDirection() << endl;
	cout << static_cast<const TImage *>(this->ProcessObject::GetInput(2))->GetDirection() << endl;
	Superclass::Update();
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::GenerateData() {
	cout << __PRETTY_FUNCTION__ << endl;
	typename TVectorImage::ConstPointer spgrData = this->GetInput();
	typename TImage::ConstPointer maskData = this->GetMask();
	typename TImage::ConstPointer B1Data   = this->GetB1();

	typename TImage::Pointer T1Data = this->GetOutputT1();
	typename TImage::Pointer PDData = this->GetOutputPD();
	typename TImage::Pointer ResData = this->GetOutputRes();

	typename TImage::RegionType region = maskData->GetLargestPossibleRegion();
	T1Data->SetRegions(region);
	PDData->SetRegions(region);
	ResData->SetRegions(region);
	T1Data->Allocate();
	PDData->Allocate();
	ResData->Allocate();

	ImageRegionConstIterator<TVectorImage> SPGRIter(spgrData, spgrData->GetLargestPossibleRegion());
	ImageRegionConstIterator<TImage> maskIter(maskData, maskData->GetLargestPossibleRegion());
	ImageRegionConstIterator<TImage> B1Iter(B1Data, B1Data->GetLargestPossibleRegion());
	ImageRegionIterator<TImage> T1Iter(T1Data, T1Data->GetLargestPossibleRegion());
	ImageRegionIterator<TImage> PDIter(PDData, PDData->GetLargestPossibleRegion());
	ImageRegionIterator<TImage> ResIter(ResData, ResData->GetLargestPossibleRegion());

	shared_ptr<SCD> model = make_shared<SCD>();
	while(!T1Iter.IsAtEnd()) {
		if (maskIter.Get()) {
			double B1 = B1Iter.Get();
			ArrayXd localAngles = m_sequence.flip() * B1;
			double T1, PD;
			VariableLengthVector<float> signalVector = SPGRIter.Get();
			Map<const ArrayXf> signalf(signalVector.GetDataPointer(), m_sequence.size());
			ArrayXd signal = signalf.cast<double>();
			VectorXd Y = signal / localAngles.sin();
			MatrixXd X(Y.rows(), 2);
			X.col(0) = signal / localAngles.tan();
			X.col(1).setOnes();
			VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
			T1 = -m_sequence.TR() / log(b[0]);
			PD = b[1] / (1. - b[0]);
			if (m_algorithm == Algos::WLLS) {
				VectorXd W(m_sequence.size());
				for (size_t n = 0; n < m_iterations; n++) {
					W = (localAngles.sin() / (1. - (exp(-m_sequence.TR()/T1)*localAngles.cos()))).square();
					b = (X.transpose() * W.asDiagonal() * X).partialPivLu().solve(X.transpose() * W.asDiagonal() * Y);
					T1 = -m_sequence.TR() / log(b[0]);
					PD = b[1] / (1. - b[0]);
				}
			} else if (m_algorithm == Algos::NLLS) {
				/*DESPOTFunctor f(spgrSequence, Pools::One, signal.cast<complex<double>>(), B1, false, false);
				NumericalDiff<DESPOTFunctor> nDiff(f);
				LevenbergMarquardt<NumericalDiff<DESPOTFunctor>> lm(nDiff);
				lm.parameters.maxfev = m_iterations;
				VectorXd p(4);
				p << PD, T1, 0., 0.; // Don't need T2 of f0 for this (yet)
				lm.lmder1(p);
				PD = p(0); T1 = p(1);*/
			}
			VectorXd pars(5); pars << PD, T1, 0., 0., B1;
			ArrayXd theory = m_sequence.signal(model, pars).abs();
			ArrayXd resids = (signal - theory);
			/*if (all_residuals) {
				ResidsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = resids.cast<float>();
			}*/
			T1Iter.Set(static_cast<float>(T1));
			PDIter.Set(static_cast<float>(PD));
			ResIter.Set(static_cast<float>(sqrt(resids.square().sum() / resids.rows()) / PD));
		} else {
			//T1Iter.Set(0);
			//PDIter.Set(0);
			//ResIter.Set(0);
		}
		++SPGRIter; ++maskIter; ++B1Iter;
		++T1Iter; ++PDIter; ++ResIter;
	}
}
} // namespace ITK


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

typedef itk::Image<float, 3> FloatImage;
typedef itk::VectorImage<float, 3> FloatVectorImage;
typedef itk::DESPOT1Filter<FloatVectorImage, FloatImage> DESPOT1;

static bool verbose = false, prompt = true, all_residuals = false;
static size_t nIterations = 4;
static string outPrefix;
static DESPOT1::Algos algo;
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

	typedef itk::ImageFileReader<FloatImage> Reader;
	typedef itk::ImageFileReader<itk::Image<float, 4>> Reader4D;
	typedef itk::ImageFileWriter<FloatImage> Writer;

	Reader::Pointer mask = Reader::New();
	Reader::Pointer B1   = Reader::New();

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'm':
				cout << "Opening mask file " << optarg << endl;
				mask->SetFileName(optarg);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'b':
				cout << "Opening B1 file: " << optarg << endl;
				B1->SetFileName(optarg);
				break;
			case 'a':
				switch (*optarg) {
					case 'l': algo = DESPOT1::Algos::LLS;  cout << "LLS algorithm selected." << endl; break;
					case 'w': algo = DESPOT1::Algos::WLLS; cout << "WLLS algorithm selected." << endl; break;
					case 'n': algo = DESPOT1::Algos::NLLS; cout << "NLLS algorithm selected." << endl; break;
					default:
						cout << "Unknown algorithm type " << optarg << endl;
						return EXIT_FAILURE;
						break;
				} break;
			case 'i':
				nIterations = atoi(optarg);
				break;
			case 'T':
				//threads.resize(atoi(optarg));
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

	//**************************************************************************
	#pragma mark Gather SPGR data
	//**************************************************************************
	cout << "Opening SPGR file: " << argv[optind] << endl;
	Reader4D::Pointer input = Reader4D::New();
	input->SetFileName(argv[optind]);
	input->Update();

	typedef ImageToVectorFilter<float> Converter;
	cout << "Creating converter" << endl;
	Converter::Pointer convert = Converter::New();
	cout << "Setting input" << endl;
	convert->SetInput(input->GetOutput());
	cout << "Updating" << endl;
	convert->Update();

	cout << "Input direction: " << endl << input->GetOutput()->GetDirection() << endl;
	cout << "Output direction: " << endl << convert->GetOutput()->GetDirection() << endl;

	//Agilent::ProcPar pp; ReadPP(spgrFile, pp);
	SPGRSimple spgrSequence(prompt);
	if (verbose) {
		cout << spgrSequence;
		cout << "Ouput prefix will be: " << outPrefix << endl;
	}

	cout << "Checking size" << endl;
	auto size = input->GetOutput()->GetLargestPossibleRegion().GetSize();
	cout << spgrSequence.size() << " " << input->GetFileName() << " " << size[3] << endl;
	if (spgrSequence.size() != size[3]) {
		throw(std::runtime_error("Specified number of flip-angles does not match number of volumes in file: " + input->GetFileName()));
	}
	auto vectorImage = convert->GetOutput();
	cout << vectorImage->GetLargestPossibleRegion() << endl;
	cout << "Creating DESPOT1 Filter" << endl;
	DESPOT1::Pointer d1 = DESPOT1::New();
	cout << "Input direction is: " << endl << convert->GetOutput()->GetDirection() << endl;
	cout << "Setting input" << endl;
	d1->SetInput(convert->GetOutput());
	cout << "Setting mask" << endl;
	mask->Update();
	d1->SetMask(mask->GetOutput());
	cout << "Setting B1" << endl;
	B1->Update();
	d1->SetB1(B1->GetOutput());
	d1->SetSequence(spgrSequence);
	d1->SetAlgorithm(algo);
	d1->SetIterations(nIterations);
	cout << "Created filter" << endl;

	/*if (verbose) {
		clock_t loopEnd = clock();
		if (voxCount > 0)
			cout << voxCount << " unmasked voxels, CPU time per voxel was "
					  << ((loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC)) << " s, ";
		cout << "finished." << endl;
	*/

	if (verbose)
		cout << "Writing results." << endl;
	outPrefix = outPrefix + "D1_";

	cout << "Creating Writers" << endl;
	Writer::Pointer T1File = Writer::New();
	Writer::Pointer PDFile = Writer::New();
	Writer::Pointer ResFile = Writer::New();

	T1File->SetFileName(outPrefix + "T1.nii");
	PDFile->SetFileName(outPrefix + "PD.nii");
	ResFile->SetFileName(outPrefix + "residual.nii");

	T1File->SetInput(d1->GetOutputT1());
	PDFile->SetInput(d1->GetOutputPD());
	ResFile->SetInput(d1->GetOutputRes());

	cout << "PrintDirections" << endl;
	d1->PrintDirections();
	cout << "Calling update" << endl;
	d1->Update();
	cout << "Updating output files" << endl;
	T1File->Update();
	PDFile->Update();
	ResFile->Update();

	/*
	if (all_residuals) {
		outHdr.intent_name = "Residuals";
		outHdr.setDim(4, spgrSequence.size());
		outFile.setHeader(outHdr);
		outFile.open(outPrefix + "residuals" + OutExt(), Nifti::Mode::Write);
		outFile.writeVolumes(ResidsVols.begin(), ResidsVols.end(), 0, spgrSequence.size());
		outFile.close();
	}*/
	cout << "All done." << endl;
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
