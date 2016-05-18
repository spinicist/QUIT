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
#include "QI/Util.h"

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

/*
 * Helper functions
 */

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

// Complex equivalent of the dot product
template<typename T> inline T cdot(const complex<T> &a, const complex<T> &b) {
    return real(a * conj(b));
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
            QI_EXCEPTION("Must have a minimum of 4 phase-cycling patterns.");
		if ((p % 2) != 0)
            QI_EXCEPTION("Number of phases must be even.");

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
            QI_EXCEPTION("Input size and number of phases do not match");
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
                Map<const ArrayXcf> allInput(inputVector.GetDataPointer(), m_phases);
                const ArrayXcd a = allInput.head(m_lines).cast<complex<double>>();
                const ArrayXcd b = allInput.tail(m_lines).cast<complex<double>>();

                complex<double> sum(0., 0.);
                for (size_t i = 0; i < m_lines; i++) {
                    for (size_t j = i + 1; j < m_lines; j++) {
                        const complex<double> di = b[i] -  a[i], dj = b[j] - a[j];
                        const complex<double> ni(di.imag(), -di.real()), nj(dj.imag(), -dj.real());

                        const double mu = cdot(a[j] - a[i], nj) / cdot(di, nj);
                        const double nu = cdot(a[i] - a[j], ni) / cdot(dj, ni);
                        const double xi = 1.0 - pow(cdot(di, dj) / (abs(di)*abs(dj)), 2.0);

                        const complex<double> cs = (a[i] + a[j] + b[i] + b[j]) / 4.0;
                        const complex<double> gs = a[i] + mu * di;

                        switch (m_Regularise) {
                        case RegEnum::None: sum += gs; break;
                        case RegEnum::Magnitude:
                            if (norm(gs) < max(max(max(norm(a[i]), norm(a[j])), norm(b[i])), norm(b[j]))) {
                                sum += gs;
                            } else {
                                sum += cs;
                            }
                            break;
                        case RegEnum::Line:
                            if ((mu > -xi) && (mu < 1 + xi) && (nu > -xi) && (nu < 1 + xi)) {
                                sum += gs;
                            } else {
                                sum += cs;
                            }
                            break;
                        }
                    }
                }
                outputIter.Set(static_cast<complex<float>>(sum / static_cast<double>(m_crossings)));
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
            QI_EXCEPTION("Must have a minimum of 4 phase-cycling patterns.");
		if ((p % 2) != 0)
            QI_EXCEPTION("Number of phases must be even.");
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
            QI_EXCEPTION("Input size and number of phases do not match");
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

                complex<double> output = 0.;
                VariableLengthVector<complex<float>> inputVector = inputIter.GetCenterPixel();
                Map<const ArrayXcf> allInput(inputVector.GetDataPointer(), m_phases);
                ArrayXcd a_center = allInput.head(m_lines).cast<complex<double>>();
                ArrayXcd b_center = allInput.tail(m_lines).cast<complex<double>>();
                for (int i = 0; i < m_lines; i++) {
                    double num = 0., den = 0.;
                    for (int p = 0; p < inputIter.Size(); ++p) {
                        const complex<double> Id = static_cast<complex<double>>(pass1Iter.GetPixel(p));
                        VariableLengthVector<complex<float>> inputVector = inputIter.GetPixel(p);
                        Map<const VectorXcf> allInput(inputVector.GetDataPointer(), m_phases);
                        const ArrayXcd a_pixel = allInput.head(m_lines).cast<complex<double>>();
                        const ArrayXcd b_pixel = allInput.tail(m_lines).cast<complex<double>>();

                        num += real(conj(b_pixel[i] - Id)*(b_pixel[i] - a_pixel[i]) + conj(b_pixel[i] - a_pixel[i])*(b_pixel[i] - Id));
                        den += real(conj(a_pixel[i] - b_pixel[i])*(a_pixel[i] - b_pixel[i]));
                    }
                    double w = num / (2. * den);
                    if (isfinite(w)) {
                        output += w*a_center[i] + (1. - w)*b_center[i];
                    }
                }
                outputIter.Set(static_cast<complex<float>>(output / static_cast<double>(m_lines)));
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

    QI::VolumeF::Pointer mask = ITK_NULLPTR;
    auto gs = itk::GSFilter::New();
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
            case 'm':
                if (verbose) cout << "Reading mask file " << optarg << endl;
                mask = QI::ReadImage(optarg);
                break;
			case 'o':
				prefix = optarg;
				if (verbose) cout << "Output prefix will be: " << prefix << endl;
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
				cout << QI::GetVersion() << endl << usage << endl;
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
		cout << "Incorrect number of arguments." << endl << QI::GetVersion() << endl << usage << endl;
		return EXIT_FAILURE;
	}
	if (verbose) cout << "Opening input file: " << argv[optind] << endl;
	string fname(argv[optind++]);

    auto inFile = itk::ImageFileReader<QI::SeriesXF>::New();
    auto reorderVolumes = QI::ReorderSeriesXF::New();
    auto reorderPhase = QI::ReorderSeriesXF::New();

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
    auto blockVector = QI::SeriesToVectorXF::New();
    blockVector->SetInput(reorderPhase->GetOutput());
    blockVector->SetBlockSize(nPhases);

    typename itk::ImageToImageFilter<QI::VectorVolumeXF, QI::VolumeXF>::Pointer process = ITK_NULLPTR;
    switch (output) {
    case OutEnum::GS: {
        gs->SetInput(blockVector->GetOutput());
        gs->SetPhases(nPhases);
        if (mask)
            gs->SetMask(mask);
        if (do_2pass) {
            auto p2 = itk::MinEnergyFilter::New();
            p2->SetInput(blockVector->GetOutput());
            p2->SetPass1(gs->GetOutput());
            if (mask)
                p2->SetMask(mask);
            process = p2;
        } else {
            process = gs;
        }
    } break;
    case OutEnum::CS: {
        auto cs = itk::UnaryFunctorImageFilter<QI::VectorVolumeXF, QI::VolumeXF, VectorMean<complex<float>>>::New();
        cs->SetInput(blockVector->GetOutput());
        process = cs;
    } break;
    case OutEnum::MagMean: {
        auto filter = itk::UnaryFunctorImageFilter<QI::VectorVolumeXF, QI::VolumeXF, VectorMagMean<float>>::New();
        filter->SetInput(blockVector->GetOutput());
        process = filter;
    } break;
    case OutEnum::RMS: {
        auto filter = itk::UnaryFunctorImageFilter<QI::VectorVolumeXF, QI::VolumeXF, VectorRMS<complex<float>>>::New();
        filter->SetInput(blockVector->GetOutput());
        process = filter;
    } break;
    case OutEnum::Max: {
        auto filter = itk::UnaryFunctorImageFilter<QI::VectorVolumeXF, QI::VolumeXF, VectorMax<complex<float>>>::New();
        filter->SetInput(blockVector->GetOutput());
        process = filter;
    } break;
    }

    auto outTiler = itk::TileImageFilter<QI::VolumeXF, QI::SeriesXF>::New();
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
    QI::SeriesXF::Pointer output = outTiler->GetOutput();
    output->DisconnectPipeline();
    output->SetDirection(inFile->GetOutput()->GetDirection());
    output->SetSpacing(inFile->GetOutput()->GetSpacing());
    
    if (prefix == "")
        prefix = QI::StripExt(fname).append("_nobands");
    string outname = prefix;
    outname.append(OutExt());
    if (verbose) cout << "Output filename: " << outname << endl;
    if (output_magnitude) {
        QI::WriteMagnitudeImage<QI::SeriesXF>(output, outname);
    } else {
        QI::WriteImage<QI::SeriesXF>(output, outname);
    }
    if (verbose) cout << "Finished." << endl;
	return EXIT_SUCCESS;
}
