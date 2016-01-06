/*
 *  qissfpbands.cpp
 *
 *  Created by Tobias Wood on 14/03/2014.
 *  Copyright (c) 2014 Tobias Wood.
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

#include "itkConstNeighborhoodIterator.h"

#include "Filters/ReorderVectorFilter.h"
#include "Util.h"

using namespace std;
using namespace Eigen;
using namespace QI;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: ssfpgs [options] input \n\
\n\
Input must be a single complex image with at least two pairs of opposing\n\
phase-cycles along the 4th dimension.\n\
\n\
Options:\n\
    --help, -h       : Print this message.\n\
    --verbose, -v    : Print more information.\n\
    --out, -o path   : Specify an output filename (default image base).\n\
    --mask, -m file  : Mask input with specified file.\n\
    --alt_order      : Opposing phase-cycles alternate (default is two blocks).\n\
    --threads, -T N  : Use N threads (default=hardware limit).\n\
    --save, -sL      : Save the line regularised GS (default)\n\
              M      : Save the magnitude regularised GS\n\
              G      : Save the unregularised GS\n\
              C      : Save the CS\n\
    --secondpass, -2 : Perform a 2nd pass as per Xiang and Hoff\n\
Options for multiple output volumes:\n\
    --ph_incs, -p N : Number of phase increments (default is 4).\n\
    --ph_order      : Data order is phase, then flip-angle (default opposite).\n"
};

bool verbose = false, do_pass2 = false;
int order_phase = false, order_alternate = false;
static size_t nPhases = 4;
static string prefix;
const struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
    {"alt_order", no_argument, &order_alternate, true},
    {"ph_order", required_argument, &order_phase, true},
    {"ph_incs", required_argument, 0, 'p'},
	{"threads", required_argument, 0, 'T'},
	{"save", required_argument, 0, 's'},
    {"secondpass", no_argument, 0, '2'},
	{0, 0, 0, 0}
};
const char *short_options = "hvo:m:s:p:T:2";
// From Knuth, surprised this isn't in STL
unsigned long long choose(unsigned long long n, unsigned long long k) {
	if (k > n)
		return 0;

	unsigned long long r = 1;
	for (unsigned long long d = 1; d <= k; ++d) {
		r *= n--;
		r /= d;
	}
	return r;
}

namespace itk {

class GSFilter : public ImageToImageFilter<VectorImage<complex<float>, 3>, Image<complex<float>, 3>>
{
public:
    enum class Save { LR, MR, GS, CS };

protected:
	size_t m_flips, m_lines, m_crossings, m_phases = 0;
	Save m_mode = Save::LR;
public:
	/** Standard class typedefs. */
    typedef VectorImage<complex<float>, 3> TIn;
    typedef Image<complex<float>, 3>       TOut;
    typedef Image<float, 3>                TMask;
    typedef GSFilter                       Self;
    typedef ImageToImageFilter<TIn, TOut>  Superclass;
    typedef SmartPointer<Self>             Pointer;

	itkNewMacro(Self); /** Method for creation through the object factory. */
    itkTypeMacro(GSFilter, ImageToImageFilter); /** Run-time type information (and related methods). */

	void SetSave(Save s) { m_mode = s; }
	Save GetSave() { return m_mode; }
	void SetPhases(const size_t p) {
		if (p < 4)
			throw(runtime_error("Must have a minimum of 4 phase-cycling patterns."));
		if ((p % 2) != 0)
			throw(runtime_error("Number of phases must be even."));

		m_phases = p;
		m_lines = m_phases / 2;
		m_crossings = choose(m_lines, 2);
        this->Modified();
	}
    void SetInput(const TIn *img) override {
        this->SetNthInput(0, const_cast<TIn*>(img));
	}
	void SetMask(const TMask *mask) { this->SetNthInput(1, const_cast<TMask*>(mask)); }
    typename TIn::ConstPointer GetInput() const { return static_cast<const TIn *>(this->ProcessObject::GetInput(0)); }
	typename TMask::ConstPointer GetMask() const { return static_cast<const TMask *>(this->ProcessObject::GetInput(1)); }

	virtual void GenerateOutputInformation() override {
		Superclass::GenerateOutputInformation();
		if ((this->GetInput()->GetNumberOfComponentsPerPixel() % m_phases) != 0) {
			throw(std::runtime_error("Input size and number of phases do not match"));
		}
		m_flips = (this->GetInput()->GetNumberOfComponentsPerPixel() / m_phases);
		auto op = this->GetOutput();
		op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
		op->SetNumberOfComponentsPerPixel(m_flips);
		op->Allocate();
	}

protected:
    GSFilter() {
		this->SetNumberOfRequiredInputs(1);
		this->SetNumberOfRequiredOutputs(1);
		this->SetNthOutput(0, this->MakeOutput(0));
		this->SetPhases(4);
	}
    ~GSFilter() {}

	DataObject::Pointer MakeOutput(unsigned int idx) {
		//std::cout <<  __PRETTY_FUNCTION__ << endl;
		if (idx == 0) {
            DataObject::Pointer output = (TOut::New()).GetPointer();
			return output.GetPointer();
		} else {
			std::cerr << "No output " << idx << std::endl;
			return NULL;
		}
	}

    virtual void ThreadedGenerateData(const TIn::RegionType &region, ThreadIdType threadId) override {
		//std::cout <<  __PRETTY_FUNCTION__ << endl;
        ImageRegionConstIterator<TIn> inputIter(this->GetInput(), region);
		auto m = this->GetMask();
		ImageRegionConstIterator<TMask> maskIter;
		if (m) {
			maskIter = ImageRegionConstIterator<TMask>(m, region);
		}
        ImageRegionIterator<TOut> outputIter(this->GetOutput(), region);

		while(!inputIter.IsAtEnd()) {
			if (!m || maskIter.Get()) {
				VariableLengthVector<complex<float>> inputVector = inputIter.Get();
				VariableLengthVector<complex<float>> outputVector(m_flips);

                Map<const ArrayXcf> allInput(inputVector.GetDataPointer(), m_phases);

                ArrayXcd a = allInput.head(m_lines).cast<complex<double>>();
                ArrayXcd b = allInput.tail(m_lines).cast<complex<double>>();

                Eigen::MatrixXd sols = Eigen::MatrixXd::Zero(2, m_crossings);
                size_t si = 0;
                for (size_t li = 0; li < m_lines; li++) {
                    for (size_t lj = li + 1; lj < m_lines; lj++) {
                        // Convert to 2D representation
                        Vector2d a_i{a(li).real(), a(li).imag()};
                        Vector2d a_j{a(lj).real(), a(lj).imag()};
                        Vector2d b_i{b(li).real(), b(li).imag()};
                        Vector2d b_j{b(lj).real(), b(lj).imag()};

                        Vector2d d_i = (b_i - a_i);
                        Vector2d d_j = (b_j - a_j);
                        Vector2d n_i{d_i[1], -d_i[0]};
                        Vector2d n_j{d_j[1], -d_j[0]};

                        double mu = (a_j - a_i).dot(n_j) / d_i.dot(n_j);
                        double nu = (a_i - a_j).dot(n_i) / d_j.dot(n_i);
                        double xi = 1.0 - pow(d_i.dot(d_j) / (d_i.norm() * d_j.norm()),2.0);

                        Vector2d cs = (a_i + a_j + b_i + b_j) / 4.0;
                        Vector2d gs = a_i + mu * d_i;

                        bool line_reg = true;
                        // Do the logic this way round so NaN does not propagate
                        if ((mu > -xi) && (mu < 1 + xi) && (nu > -xi) && (nu < 1 + xi))
                            line_reg = false;

                        bool mag_reg = true;
                        double maxnorm = max(max(max(a_i.norm(), a_j.norm()), b_i.norm()), b_j.norm());
                        if (gs.norm() < maxnorm) {
                            mag_reg = false;
                        }
                        switch (m_mode) {
                            case Save::LR: sols.col(si) = line_reg ? cs : gs; break;
                            case Save::MR: sols.col(si) = mag_reg ? cs : gs; break;
                            case Save::GS: sols.col(si) = gs; break;
                            case Save::CS: sols.col(si) = cs; break;
                        }
                        si++;
                    }
                }
                Vector2d mean_sol = sols.rowwise().mean();
                // Convert back to complex
                outputIter.Set(complex<float>(static_cast<float>(mean_sol[0]), static_cast<float>(mean_sol[1])));
			}
			++inputIter;
			if (m)
				++maskIter;
			++outputIter;
		}
	}

private:
    GSFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
};

class TwoPassGSFilter : public ImageToImageFilter<VectorImage<complex<float>, 3>, Image<complex<float>, 3>>
{
protected:
	size_t m_flips, m_phases, m_lines = 0;

public:
    typedef VectorImage<complex<float>, 3>     TInputImage;
    typedef Image<complex<float>, 3>           TOutputImage;
	typedef Image<float, 3>                    TMask;
    typedef TwoPassGSFilter                   Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
	typedef SmartPointer<Self>                 Pointer;

	itkNewMacro(Self);
    itkTypeMacro(TwoPassGSFilter, ImageToImageFilter);

	void SetPhases(const size_t p) {
		if (p < 4)
			throw(runtime_error("Must have a minimum of 4 phase-cycling patterns."));
		if ((p % 2) != 0)
			throw(runtime_error("Number of phases must be even."));
		m_phases = p;
		m_lines = m_phases / 2;
        this->Modified();
	}
    void SetInput(const TInputImage *img) override { this->SetNthInput(0, const_cast<TInputImage*>(img)); }
    void SetPass1(const TOutputImage *img) { this->SetNthInput(1, const_cast<TOutputImage*>(img)); }
	void SetMask(const TMask *mask) { this->SetNthInput(2, const_cast<TMask*>(mask)); }
    typename TInputImage::ConstPointer GetInput() const { return static_cast<const TInputImage *>(this->ProcessObject::GetInput(0)); }
    typename TOutputImage::ConstPointer GetPass1() const { return static_cast<const TOutputImage *>(this->ProcessObject::GetInput(1)); }
	typename TMask::ConstPointer GetMask() const { return static_cast<const TMask *>(this->ProcessObject::GetInput(2)); }

	virtual void GenerateOutputInformation() override {
		Superclass::GenerateOutputInformation();
		if ((this->GetInput()->GetNumberOfComponentsPerPixel() % m_phases) != 0) {
			throw(std::runtime_error("Input size and number of phases do not match"));
		}
		m_flips = (this->GetInput()->GetNumberOfComponentsPerPixel() / m_phases);
		auto op = this->GetOutput();
		op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
		op->SetNumberOfComponentsPerPixel(m_flips);
		op->Allocate();
	}

protected:
    TwoPassGSFilter() {
		this->SetNumberOfRequiredInputs(2);
		this->SetNumberOfRequiredOutputs(1);
		this->SetNthOutput(0, this->MakeOutput(0));
		this->SetPhases(4);
	}
    ~TwoPassGSFilter() {}

	DataObject::Pointer MakeOutput(unsigned int idx) {
		//std::cout <<  __PRETTY_FUNCTION__ << endl;
		if (idx == 0) {
            DataObject::Pointer output = (TOutputImage::New()).GetPointer();
			return output.GetPointer();
		} else {
			std::cerr << "No output " << idx << std::endl;
			return NULL;
		}
	}

    virtual void ThreadedGenerateData(const TInputImage::RegionType &region, ThreadIdType threadId) override {
		//std::cout <<  __PRETTY_FUNCTION__ << endl;
        ConstNeighborhoodIterator<TInputImage>::RadiusType radius;
		radius.Fill(1);
        ConstNeighborhoodIterator<TInputImage> inputIter(radius, this->GetInput(), region);
        ConstNeighborhoodIterator<TOutputImage> pass1Iter(radius, this->GetPass1(), region);

		auto m = this->GetMask();
		ImageRegionConstIterator<TMask> maskIter;
		if (m) {
			maskIter = ImageRegionConstIterator<TMask>(m, region);
		}
        ImageRegionIterator<TOutputImage> outputIter(this->GetOutput(), region);
		while(!inputIter.IsAtEnd()) {
			if (!m || maskIter.Get()) {
                double num = 0., den = 0.;

				for (int p = 0; p < inputIter.Size(); ++p) {
                    VariableLengthVector<complex<float>> inputVector = inputIter.GetPixel(p);
                    Map<const VectorXcf> allInput(inputVector.GetDataPointer(), m_phases);
                    ArrayXcd a = allInput.head(m_lines).cast<complex<double>>();
                    ArrayXcd b = allInput.tail(m_lines).cast<complex<double>>();
                    complex<double> s_i = static_cast<complex<double>>(pass1Iter.GetPixel(p));
                    for (int li = 0; li < m_lines; li++) {
                        complex<double> a_i = a[li];
                        complex<double> b_i = b[li];
                        num += real(conj(b_i - s_i)*(a_i - b_i) + conj(a_i - b_i)*(b_i - s_i));
                        den += real(conj(a_i - b_i)*(a_i - b_i));
                    }
                }

                double w = -num / (2. * den);
				VariableLengthVector<complex<float>> inputVector = inputIter.GetCenterPixel();
                Map<const ArrayXcf> allInput(inputVector.GetDataPointer(), m_phases);
                complex<float> output = 0.;
                if (isfinite(w)) {
                    ArrayXcd a = allInput.head(m_lines).cast<complex<double>>();
                    ArrayXcd b = allInput.tail(m_lines).cast<complex<double>>();
                    for (int li = 0; li < m_lines; li++) {
                        complex<double> s = (a[li]*w + (1.f-w)*b[li]);
                        output += static_cast<complex<float>>(s / static_cast<double>(m_lines));
                    }
                }
                outputIter.Set(output);
			}
			++inputIter;
			++pass1Iter;
			if (m)
				++maskIter;
			++outputIter;
		}
		//std::cout << "End " << __PRETTY_FUNCTION__ << std::endl;
	}

private:
    TwoPassGSFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
};

} // End namespace itk

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();

	itk::ImageFileReader<itk::Image<float, 3>>::Pointer mask = ITK_NULLPTR;
    auto pass1 = itk::GSFilter::New();
    auto pass2 = itk::TwoPassGSFilter::New();
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'm':
				if (verbose) cout << "Reading mask file " << optarg << endl;
				mask = itk::ImageFileReader<itk::Image<float, 3>>::New();
				mask->SetFileName(optarg);
				pass1->SetMask(mask->GetOutput());
				pass2->SetMask(mask->GetOutput());
				break;
			case 'o':
				prefix = optarg;
				cout << "Output prefix will be: " << prefix << endl;
				break;
            case 'F': order_phase = true; break;
			case 'p': nPhases = atoi(optarg); break;
            case 'a': order_alternate = true; break;
			case 's':
				switch(*optarg) {
                    case 'L': pass1->SetSave(itk::GSFilter::Save::LR); break;
                    case 'M': pass1->SetSave(itk::GSFilter::Save::MR); break;
                    case 'G': pass1->SetSave(itk::GSFilter::Save::GS); break;
                    case 'C': pass1->SetSave(itk::GSFilter::Save::CS); break;
					default:
						cerr << "Unknown desired save image: " << *optarg << endl;
						return EXIT_FAILURE;
				}
				break;
            case '2': do_pass2 = true; break;
			case 'T':
				itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg));
				break;
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

    auto inFile = QI::ReadTimeseriesXF::New();
    auto reorderVolumes = QI::ReorderTimeseriesXF::New();
    auto reorderPhase = QI::ReorderTimeseriesXF::New();

    inFile->SetFileName(fname);
    inFile->Update(); // We need to know the number of input volumes to work out the number of output volumes
    size_t nVols = inFile->GetOutput()->GetLargestPossibleRegion().GetSize()[3] / nPhases;
    if (verbose) {
        cout << "Number of phase increments is " << nPhases << endl;
        cout << "Number of volumes to process is " << nVols << endl;
    }
    // Re-order once to get all phase-increments for one output volume together
    reorderVolumes->SetInput(inFile->GetOutput());
    if (!order_phase) {
        reorderVolumes->SetStride(nVols);
    }
    // Re-order again within blocks to separate opposing phase-increments
    reorderPhase->SetInput(reorderVolumes->GetOutput());
    if (order_alternate) {
        reorderPhase->SetStride(2);
        reorderPhase->SetBlockSize(nPhases);
    }
    reorderPhase->Update();

    if (verbose) cout << "Reordered data" << endl;
    auto blockVector = QI::TimeseriesToVectorXF::New();
    blockVector->SetInput(reorderPhase->GetOutput());
    blockVector->SetBlockSize(nPhases);

    if (verbose) cout << "Setting up passes" << endl;
    pass1->SetPhases(nPhases);
    pass1->SetInput(blockVector->GetOutput());
    pass2->SetPhases(nPhases);
    pass2->SetInput(blockVector->GetOutput());
    pass2->SetPass1(pass1->GetOutput());

    auto pass1Tiler = itk::TileImageFilter<QI::ImageXF, QI::TimeseriesXF>::New();
    auto pass2Tiler = itk::TileImageFilter<QI::ImageXF, QI::TimeseriesXF>::New();
    itk::FixedArray<unsigned int, 4> layout;
    layout[0] = layout[1] = layout[2] = 1; layout[3] = nVols;
    pass1Tiler->SetLayout(layout);
    pass2Tiler->SetLayout(layout);

    for (size_t i = 0; i < nVols; i++) {
        if (verbose) cout << "Processing output volume " << i << endl;
        blockVector->SetBlockStart(i * nPhases);
        blockVector->Update();

        pass1->Update();
        auto pass1Image = pass1->GetOutput();
        pass1Tiler->SetInput(i, pass1Image);
        pass1Image->DisconnectPipeline();

        if (do_pass2) {
            pass2->Update();
            auto pass2Image = pass2->GetOutput();
            pass2Tiler->SetInput(i, pass2Image);
            pass2Image->DisconnectPipeline();
        }
    }
    pass1Tiler->UpdateLargestPossibleRegion();
    if (do_pass2)
        pass2Tiler->UpdateLargestPossibleRegion();

	auto outFile = QI::WriteTimeseriesXF::New();
	if (prefix == "")
        prefix = QI::StripExt(fname);
	string outname = prefix;
	switch (pass1->GetSave()) {
        case itk::GSFilter::Save::LR: outname.append("_lreg"); break;
        case itk::GSFilter::Save::MR: outname.append("_mreg"); break;
        case itk::GSFilter::Save::GS: outname.append("_gs"); break;
        case itk::GSFilter::Save::CS: outname.append("_cs"); break;
	}
	if (do_pass2) {
		outname.append("_2p");
        outFile->SetInput(pass2Tiler->GetOutput());
	} else {
        outFile->SetInput(pass1Tiler->GetOutput());
	}
	outname.append(OutExt());
	if (verbose) cout << "Output filename: " << outname << endl;
	outFile->SetFileName(outname);
	if (verbose) cout << "Processing..." << endl;
	outFile->Update();
	if (verbose) cout << "Finished." << endl;
	return EXIT_SUCCESS;
}
