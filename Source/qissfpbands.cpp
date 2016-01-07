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

#include "itkUnaryFunctorImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
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
    --magnitude, -M  : Output a magnitude image (default is complex).\n\
Output options (mutually exclusive):\n\
    --magmean        : Output the mean of the magnitudes of the phase increments.\n\
    --rms            : Output the root mean square.\n\
    --max            : Output the maximum intensity projection.\n\
    --cs             : Output the Complex Solution (complex mean).\n\
    --gs             : Output the Geometric Solution (default).\n\
Regularisation options for Geometric Solution\n\
    --regularise, -R M : Use magnitude regularisation (from Xiang and Hoff).\n\
                     L : Use line regularisation (default).\n\
                     N : Don't regularise (none).\n\
    --secondpass, -2   : Use the energy-minimisation scheme from Xiang and Hoff.\n\
Options for multiple output volumes:\n\
    --ph_incs, -p N : Number of phase increments (default is 4).\n\
    --ph_order      : Data order is phase, then flip-angle (default opposite).\n"
};

enum OutEnum { MagMean = 1, RMS, Max, CS, GS };
bool verbose = false, output_magnitude = false;
int order_phase = false, order_alternate = false, output = OutEnum::GS, do_2pass;
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
    {"magnitude", no_argument, 0, 'M'},
    {"magmean", no_argument, &output, OutEnum::MagMean},
    {"rms", no_argument, &output, OutEnum::RMS},
    {"max", no_argument, &output, OutEnum::Max},
    {"cs", no_argument, &output, OutEnum::CS},
    {"gs", no_argument, &output, OutEnum::GS},
    {"regularise", required_argument, 0, 'R'},
    {"secondpass", no_argument, &do_2pass, true},
	{0, 0, 0, 0}
};
const char *short_options = "hvo:m:s:p:T:MR:2";
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

class GSFilter : public ImageToImageFilter<VectorImage<complex<float>, 3>, Image<complex<float>, 3>> {
public:
    enum class RegEnum { None = 0, Line, Magnitude };
protected:
	size_t m_flips, m_lines, m_crossings, m_phases = 0;
    RegEnum m_Regularise = RegEnum::Line;
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

    itkSetMacro(Regularise, RegEnum);
    itkGetMacro(Regularise, RegEnum);

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
    void SetInput(const TIn *img) override { this->SetNthInput(0, const_cast<TIn*>(img)); }
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
                        switch (m_Regularise) {
                            case RegEnum::Line:      sols.col(si) = line_reg ? cs : gs; break;
                            case RegEnum::Magnitude: sols.col(si) = mag_reg ? cs : gs; break;
                            case RegEnum::None:      sols.col(si) = gs; break;
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

class MinEnergyFilter : public ImageToImageFilter<VectorImage<complex<float>, 3>, Image<complex<float>, 3>>
{
protected:
	size_t m_flips, m_phases, m_lines = 0;

public:
    typedef VectorImage<complex<float>, 3>     TInputImage;
    typedef Image<complex<float>, 3>           TOutputImage;
	typedef Image<float, 3>                    TMask;
    typedef MinEnergyFilter                   Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
	typedef SmartPointer<Self>                 Pointer;

	itkNewMacro(Self);
    itkTypeMacro(MinEnergyFilter, ImageToImageFilter);

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
    MinEnergyFilter() {
		this->SetNumberOfRequiredInputs(2);
		this->SetNumberOfRequiredOutputs(1);
		this->SetNthOutput(0, this->MakeOutput(0));
		this->SetPhases(4);
	}
    ~MinEnergyFilter() {}

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
    MinEnergyFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
};

} // End namespace itk

template<class T> class VectorMean {
public:
    VectorMean() {};
    ~VectorMean() {};
    bool operator!=( const VectorMean & ) const { return false; }
    bool operator==( const VectorMean &other ) const { return !(*this != other); }
    inline T operator()( const itk::VariableLengthVector<T> & v ) const {
        T sum = T();
        for (size_t i = 0; i < v.Size(); i++) {
            sum += v[i];
        }
        return sum / static_cast<T>(v.Size());
    }
};

template<class T> class VectorMagMean {
public:
    VectorMagMean() {};
    ~VectorMagMean() {};
    bool operator!=( const VectorMagMean & ) const { return false; }
    bool operator==( const VectorMagMean &other ) const { return !(*this != other); }
    inline complex<T> operator()( const itk::VariableLengthVector<complex<T>> & v ) const {
        T sum = T();
        for (size_t i = 0; i < v.Size(); i++) {
            sum += std::abs(v[i]);
        }
        return complex<T>(sum / v.Size(), 0.);
    }
};

template<class T> class VectorRMS {
public:
    VectorRMS() {};
    ~VectorRMS() {};
    bool operator!=( const VectorRMS & ) const { return false; }
    bool operator==( const VectorRMS &other ) const { return !(*this != other); }
    inline T operator()( const itk::VariableLengthVector<T> & v ) const {
        T sum = T();
        for (size_t i = 0; i < v.Size(); i++) {
            sum += v[i] * v[i];
        }
        return std::sqrt(sum / T(v.Size()));
    }
};

template<class T> class VectorMax {
public:
    VectorMax() {};
    ~VectorMax() {};
    bool operator!=( const VectorMax & ) const { return false; }
    bool operator==( const VectorMax &other ) const { return !(*this != other); }
    inline T operator()( const itk::VariableLengthVector<T> & v ) const {
        T max = std::numeric_limits<T>::lowest();
        for (size_t i = 0; i < v.Size(); i++) {
            if (std::abs(v[i]) > std::abs(max)) max = v[i];
        }
        return max;
    }
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();

    QI::ReadImageF::Pointer mask = ITK_NULLPTR;
    auto gs = itk::GSFilter::New();
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
            case 'F': order_phase = true; break;
			case 'p': nPhases = atoi(optarg); break;
            case 'a': order_alternate = true; break;
            case 'M': output_magnitude = true; break;
            case 'R':
                switch(*optarg) {
                    case 'L': gs->SetRegularise(itk::GSFilter::RegEnum::Line); break;
                    case 'M': gs->SetRegularise(itk::GSFilter::RegEnum::Magnitude); break;
                    case 'N': gs->SetRegularise(itk::GSFilter::RegEnum::None); break;
                    default:
                        cerr << "Unknown regularisation strategy: " << *optarg << endl;
                        return EXIT_FAILURE;
				}
				break;
            case '2': do_2pass = true; break;
			case 'T':
				itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg));
				break;
			case 'h':
				cout << usage << endl;
				return EXIT_SUCCESS;
            case 0:   // getopt set a flag from the long_opts
                break;
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

    typename itk::ImageToImageFilter<QI::VectorImageXF, QI::ImageXF>::Pointer process = ITK_NULLPTR;
    switch (output) {
    case OutEnum::GS: {
        gs->SetInput(blockVector->GetOutput());
        gs->SetPhases(nPhases);
        if (mask)
            gs->SetMask(mask->GetOutput());
        if (do_2pass) {
            auto p2 = itk::MinEnergyFilter::New();
            p2->SetInput(blockVector->GetOutput());
            p2->SetPass1(gs->GetOutput());
            if (mask)
                p2->SetMask(mask->GetOutput());
            process = p2;
        } else {
            process = gs;
        }
    } break;
    case OutEnum::CS: {
        auto cs = itk::UnaryFunctorImageFilter<QI::VectorImageXF, QI::ImageXF, VectorMean<complex<float>>>::New();
        cs->SetInput(blockVector->GetOutput());
        process = cs;
    } break;
    case OutEnum::MagMean: {
        auto filter = itk::UnaryFunctorImageFilter<QI::VectorImageXF, QI::ImageXF, VectorMagMean<float>>::New();
        filter->SetInput(blockVector->GetOutput());
        process = filter;
    } break;
    case OutEnum::RMS: {
        auto filter = itk::UnaryFunctorImageFilter<QI::VectorImageXF, QI::ImageXF, VectorRMS<complex<float>>>::New();
        filter->SetInput(blockVector->GetOutput());
        process = filter;
    } break;
    case OutEnum::Max: {
        auto filter = itk::UnaryFunctorImageFilter<QI::VectorImageXF, QI::ImageXF, VectorMax<complex<float>>>::New();
        filter->SetInput(blockVector->GetOutput());
        process = filter;
    } break;
    }

    auto outTiler = itk::TileImageFilter<QI::ImageXF, QI::TimeseriesXF>::New();
    itk::FixedArray<unsigned int, 4> layout;
    layout[0] = layout[1] = layout[2] = 1; layout[3] = nVols;
    outTiler->SetLayout(layout);

    for (size_t i = 0; i < nVols; i++) {
        if (verbose) cout << "Processing output volume " << i << endl;
        blockVector->SetBlockStart(i * nPhases);
        //blockVector->Update();
        process->Update();
        auto volume = process->GetOutput();
        outTiler->SetInput(i, volume);
        volume->DisconnectPipeline();

    }
    outTiler->UpdateLargestPossibleRegion();

    if (prefix == "")
        prefix = QI::StripExt(fname).append("_nobands");
    string outname = prefix;
    outname.append(OutExt());
    if (verbose) cout << "Output filename: " << outname << endl;
    if (output_magnitude) {
        auto mag = itk::ComplexToModulusImageFilter<QI::TimeseriesXF, QI::TimeseriesF>::New();
        auto outFile = QI::WriteTimeseriesF::New();
        mag->SetInput(outTiler->GetOutput());
        outFile->SetInput(mag->GetOutput());
        outFile->SetFileName(outname);
        outFile->Update();
    } else {
        auto outFile = QI::WriteTimeseriesXF::New();
        outFile->SetInput(outTiler->GetOutput());
        outFile->SetFileName(outname);
        outFile->Update();
    }
    if (verbose) cout << "Finished." << endl;
	return EXIT_SUCCESS;
}
