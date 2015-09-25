/*
 *  qiunwrap.cpp
 *
 *  Created by Tobias Wood on 11/06/2015.
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
#include "Eigen/Dense"

#include "itkImageSource.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMultiplyImageFilter.h"
#include "itkForwardFFTImageFilter.h"
#include "itkInverseFFTImageFilter.h"
#include "itkFFTPadImageFilter.h"
#include "Types.h"
#include "Util.h"

using namespace std;
using namespace Eigen;

namespace itk {

class DiscreteLaplacePhaseFilter : public ImageToImageFilter<QI::ImageF, QI::ImageF> {
protected:

public:
	/** Standard class typedefs. */
    typedef QI::ImageF     TImage;

	typedef DiscreteLaplacePhaseFilter         Self;
	typedef ImageToImageFilter<TImage, TImage> Superclass;
	typedef SmartPointer<Self>                 Pointer;
	typedef typename TImage::RegionType        RegionType;

	itkNewMacro(Self);
	itkTypeMacro(DiscreteLaplacePhaseFilter, DiscreteLaplacePhaseFilter);

    void SetInput(const TImage *img) override      { this->SetNthInput(0, const_cast<TImage*>(img)); }
    void SetMask(const TImage *img)                { this->SetNthInput(1, const_cast<TImage*>(img)); }
	typename TImage::ConstPointer GetInput() const { return static_cast<const TImage *>(this->ProcessObject::GetInput(0)); }
	typename TImage::ConstPointer GetMask() const  { return static_cast<const TImage *>(this->ProcessObject::GetInput(1)); }

	virtual void GenerateOutputInformation() override {
		Superclass::GenerateOutputInformation();
		auto op = this->GetOutput();
		op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
		op->Allocate();
	}

protected:
	DiscreteLaplacePhaseFilter() {
		this->SetNumberOfRequiredInputs(1);
		this->SetNumberOfRequiredOutputs(1);
		this->SetNthOutput(0, this->MakeOutput(0));
	}
	~DiscreteLaplacePhaseFilter() {}

	DataObject::Pointer MakeOutput(unsigned int idx) {
		//std::cout <<  __PRETTY_FUNCTION__ << endl;
		if (idx == 0) {
			DataObject::Pointer output = (TImage::New()).GetPointer();
			return output.GetPointer();
		} else {
			std::cerr << "No output " << idx << std::endl;
			return NULL;
		}
	}

    virtual void ThreadedGenerateData(const RegionType &region, ThreadIdType threadId) override {
		//std::cout <<  __PRETTY_FUNCTION__ << endl;
		ConstNeighborhoodIterator<TImage>::RadiusType radius;
		radius.Fill(1);
		ConstNeighborhoodIterator<TImage> inputIter(radius, this->GetInput(), region);
		ImageRegionIterator<TImage> outputIter(this->GetOutput(), region);

		ImageRegionConstIterator<TImage> maskIter;
		const auto mask = this->GetMask();
		if (mask) {
			maskIter = ImageRegionConstIterator<TImage>(mask, region);
		}

		vector<ConstNeighborhoodIterator<TImage>::OffsetType> back, fwrd;
		back.push_back({{-1, 0, 0}});
		fwrd.push_back({{ 1, 0, 0}});
		back.push_back({{ 0,-1, 0}});
		fwrd.push_back({{ 0, 1, 0}});
		back.push_back({{ 0, 0,-1}});
		fwrd.push_back({{ 0, 0, 1}});

		TImage::SpacingType spacing = this->GetInput()->GetSpacing();
		TImage::SpacingType s_sqr = spacing * spacing;
		while(!inputIter.IsAtEnd()) {
			double sum = 0;
			if (!mask || maskIter.Get()) {
				double cphase = inputIter.GetCenterPixel();
				complex<double> c = std::polar(1., cphase);
				for (int i = 0; i < fwrd.size(); ++i) {
					double bphase = inputIter.GetPixel(back[i]);
					double fphase = inputIter.GetPixel(fwrd[i]);
					complex<double> b = std::polar(1., bphase);
					complex<double> f = std::polar(1., fphase);
					sum += std::arg((f*b)/(c*c)) / s_sqr[i];
				}
			}
			outputIter.Set(sum / 7.);
			++inputIter;
			++outputIter;
			if (mask)
				++maskIter;
		}
	}

private:
	DiscreteLaplacePhaseFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
};

class DiscreteInverseLaplace : public ImageSource<QI::ImageF> {
public:
    typedef QI::ImageF             TImage;
	typedef DiscreteInverseLaplace Self;
	typedef ImageSource<TImage>    Superclass;
	typedef SmartPointer<Self>     Pointer;

	itkNewMacro(Self);
	itkTypeMacro(DiscreteInverseLaplace, ImageSource);

    void SetImageProperties(const TImage *img) {
        m_region = img->GetLargestPossibleRegion();
		m_spacing = img->GetSpacing();
		m_direction = img->GetDirection();
		m_origin = img->GetOrigin();
	}

protected:
    typename TImage::RegionType    m_region;
	typename TImage::SpacingType   m_spacing;
	typename TImage::DirectionType m_direction;
	typename TImage::PointType     m_origin;

	DiscreteInverseLaplace(){}
	~DiscreteInverseLaplace(){}
    virtual void GenerateData() override {
		typename TImage::Pointer output = this->GetOutput();
        output->SetRegions(m_region);
		output->Allocate();
		output->SetSpacing(m_spacing);
		output->SetDirection(m_direction);
		output->SetOrigin(m_origin);
		itk::ImageRegionIteratorWithIndex<TImage> imageIt(output,output->GetLargestPossibleRegion());
		imageIt.GoToBegin();
		imageIt.Set(0.); // There is a pole here
		++imageIt;
		while(!imageIt.IsAtEnd()) {
            auto index = imageIt.GetIndex() - m_region.GetIndex(); // Might be padded to a negative start
			double val = 0;
			for (int i = 0; i < 3; i++) {
                val += 2. - 2. * cos(index[i] * 2. * M_PI / m_region.GetSize()[i]);
			}
			val /= 7.;
			imageIt.Set(1./val);
			++imageIt;
		}
	}

private:
	DiscreteInverseLaplace(const Self &);
	void operator=(const Self &);
};

} // End namespace itk

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: qiunwrap [options] input \n\
\n\
Input is a single wrapped phase volume\n\
\n\
Options:\n\
	--help, -h        : Print this message.\n\
	--verbose, -v     : Print more information.\n\
	--out, -o path    : Specify an output filename (default image base).\n\
	--mask, -m file   : Mask input with specified file.\n\
	--threads, -T N   : Use N threads (default=hardware limit).\n\
	--lop, -l C       : Use Continuous Laplacian operators.\n\
	          D         Use Discrete Laplacian operators (default).\n"
};

const struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"threads", required_argument, 0, 'T'},
	{"lop", required_argument, 0, 'l'},
	{0, 0, 0, 0}
};
const char *short_options = "hvo:m:T:l:";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();

	bool verbose = false;
	string prefix;
    QI::ReadImageF::Pointer mask = ITK_NULLPTR;
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'm':
				if (verbose) cout << "Reading mask file " << optarg << endl;
                mask = QI::ReadImageF::New();
				mask->SetFileName(optarg);
				break;
			case 'o':
				prefix = optarg;
				cout << "Output prefix will be: " << prefix << endl;
				break;
			case 'l':
				break;
			case 'T': itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg)); break;
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
	if ((argc - optind) != 1) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		return EXIT_FAILURE;
	}
	if (verbose) cout << "Opening input file: " << argv[optind] << endl;
	string fname(argv[optind++]);
	if (prefix == "")
		prefix = fname.substr(0, fname.find(".nii"));
	string outname = prefix + "_unwrap" + QI::OutExt();
	if (verbose) cout << "Output filename: " << outname << endl;

    auto inFile = QI::ReadImageF::New();
	auto calcLaplace = itk::DiscreteLaplacePhaseFilter::New();
	inFile->SetFileName(fname);
	inFile->Update(); // Need the size info
	calcLaplace->SetInput(inFile->GetOutput());
	if (mask)
        calcLaplace->SetMask(mask->GetOutput());

    if (verbose) cout << "Padding image to valid FFT size." << endl;
    typedef itk::FFTPadImageFilter<QI::ImageF> PadFFTType;
    auto padFFT = PadFFTType::New();
    padFFT->SetInput(calcLaplace->GetOutput());
    padFFT->Update();

    if (verbose) {
        cout << "Padded image size: " << padFFT->GetOutput()->GetLargestPossibleRegion().GetSize() << endl;
        cout << "Calculating Forward FFT." << endl;
    }
    typedef itk::ForwardFFTImageFilter<QI::ImageF> FFFTType;
	auto forwardFFT = FFFTType::New();
    forwardFFT->SetInput(padFFT->GetOutput());
    forwardFFT->Update();

	if (verbose) cout << "Generating Inverse Laplace Kernel." << endl;
	auto inverseLaplace = itk::DiscreteInverseLaplace::New();
    inverseLaplace->SetImageProperties(padFFT->GetOutput());
    inverseLaplace->Update();

    if (verbose) cout << "Multiplying." << endl;
    auto mult = itk::MultiplyImageFilter<QI::ImageXF, QI::ImageF, QI::ImageXF>::New();
	mult->SetInput1(forwardFFT->GetOutput());
	mult->SetInput2(inverseLaplace->GetOutput());

	if (verbose) cout << "Inverse FFT." << endl;
    auto inverseFFT = itk::InverseFFTImageFilter<QI::ImageXF, QI::ImageF>::New();
	inverseFFT->SetInput(mult->GetOutput());
    inverseFFT->Update();

    if (verbose) cout << "Extracting original size image" << endl;
    auto extract = itk::ExtractImageFilter<QI::ImageF, QI::ImageF>::New();
    extract->SetInput(inverseFFT->GetOutput());
    extract->SetDirectionCollapseToSubmatrix();
    extract->SetExtractionRegion(calcLaplace->GetOutput()->GetLargestPossibleRegion());
    extract->Update();

    auto outFile = QI::WriteImageF::New();
    outFile->SetInput(extract->GetOutput());
	outFile->SetFileName(outname);
    if (verbose) cout << "Writing output." << endl;
	outFile->Update();
	if (verbose) cout << "Finished." << endl;
	return EXIT_SUCCESS;
}
