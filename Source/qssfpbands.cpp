/*
 *  ssfpbands_main.cpp
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
	--help, -h        : Print this message.\n\
	--verbose, -v     : Print more information.\n\
	--out, -o path    : Specify an output filename (default image base).\n\
	--mask, -m file   : Mask input with specified file.\n\
	--phaseflip, -F   : Data order is phase, then flip-angle (default opposite).\n\
	--phases, -p N    : Number of phase-cycling patterns used (default is 4).\n\
	--alternate, -a   : Opposing phase-cycles alternate (default is two blocks).\n\
	--threads, -T N   : Use N threads (default=hardware limit).\n\
	--save, -sL       : Save the line regularised GS (default)\n\
	          M       : Save the magnitude regularised GS\n\
	          G       : Save the unregularised GS\n\
	          C       : Save the CS\n\
	          P       : Save the PS\n\
	--secondpass, -2  : Perform a 2nd pass as per Xiang and Hoff\n"
};

static bool verbose = false, do_pass2 = false, phaseflip = false, alternate = false;
static size_t nPhases = 4;
static string prefix;
const struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"phaseflip", required_argument, 0, 'F'},
	{"phases", required_argument, 0, 'p'},
	{"alternate", no_argument, 0, 'a'},
	{"threads", required_argument, 0, 'T'},
	{"save", required_argument, 0, 's'},
	{"secondpass", optional_argument, 0, '2'},
	{0, 0, 0, 0}
};
const char *short_options = "hvo:m:Fs:p:aT:2::";
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

class FirstPassFilter : public ImageToImageFilter<VectorImage<complex<float>, 3>, VectorImage<complex<float>, 3>>
{
public:
	enum class Save { LR, MR, GS, CS, PS };

protected:
	size_t m_flips, m_lines, m_crossings, m_phases = 0;
	Save m_mode = Save::LR;
public:
	/** Standard class typedefs. */
	typedef VectorImage<complex<float>, 3>     TImage;
	typedef Image<float, 3>                    TMask;
	typedef FirstPassFilter                    Self;
	typedef ImageToImageFilter<TImage, TImage> Superclass;
	typedef SmartPointer<Self>                 Pointer;
	typedef typename TImage::RegionType        RegionType;

	itkNewMacro(Self); /** Method for creation through the object factory. */
	itkTypeMacro(FirstPassFilter, ImageToImageFilter); /** Run-time type information (and related methods). */

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
	}
	void SetInput(const TImage *img) {
		this->SetNthInput(0, const_cast<TImage*>(img));
	}
	void SetMask(const TMask *mask) { this->SetNthInput(1, const_cast<TMask*>(mask)); }
	typename TImage::ConstPointer GetInput() const { return static_cast<const TImage *>(this->ProcessObject::GetInput(0)); }
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
	FirstPassFilter() {
		this->SetNumberOfRequiredInputs(1);
		this->SetNumberOfRequiredOutputs(1);
		this->SetNthOutput(0, this->MakeOutput(0));
		this->SetPhases(4);
	}
	~FirstPassFilter() {}

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

	virtual void ThreadedGenerateData(const RegionType &region, ThreadIdType threadId) {
		//std::cout <<  __PRETTY_FUNCTION__ << endl;
		ImageRegionConstIterator<TImage> inputIter(this->GetInput(), region);
		auto m = this->GetMask();
		ImageRegionConstIterator<TMask> maskIter;
		if (m) {
			maskIter = ImageRegionConstIterator<TMask>(m, region);
		}
		ImageRegionIterator<TImage> outputIter(this->GetOutput(), region);

		while(!inputIter.IsAtEnd()) {
			if (!m || maskIter.Get()) {
				VariableLengthVector<complex<float>> inputVector = inputIter.Get();
				VariableLengthVector<complex<float>> outputVector(m_flips);

				Map<const ArrayXXcf> allInput(inputVector.GetDataPointer(), m_phases, m_flips);
				for (int f = m_flips - 1; f > -1; f--) {
					ArrayXcd a = allInput.col(f).segment(0, m_lines).cast<complex<double>>();
					ArrayXcd b = allInput.col(f).segment(m_lines, m_lines).cast<complex<double>>();

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

							Vector2d ps;
							if (f < (m_flips - 1)) { // Use the phase of the last flip-angle for regularisation
								double phase = arg(outputVector[m_flips-1]);
								Vector2d d_p{cos(phase),sin(phase)};
								double lm_i = (a_i).dot(n_i) / d_p.dot(n_i);
								double lm_j = (a_j).dot(n_j) / d_p.dot(n_j);
								Vector2d p_i = lm_i * d_p;
								Vector2d p_j = lm_j * d_p;
								ps = (p_i + p_j) / 2.0;
							}

							bool line_reg = true;
							// Do the logic this way round so NaN does not propagate
							if ((mu > -xi) && (mu < 1 + xi) && (nu > -xi) && (nu < 1 + xi))
								line_reg = false;

							bool mag_reg = true;
							double maxnorm = max(max(max(a_i.norm(), a_j.norm()), b_i.norm()), b_j.norm());
							if (gs.norm() < maxnorm) {
								mag_reg = false;
							}
							if (ps.norm() > maxnorm) {
								ps = cs;
							}
							switch (m_mode) {
								case Save::LR: sols.col(si) = line_reg ? ps : gs; break;
								case Save::MR: sols.col(si) = mag_reg ? ps : gs; break;
								case Save::GS: sols.col(si) = gs; break;
								case Save::PS: sols.col(si) = ps; break;
								case Save::CS: sols.col(si) = cs; break;
							}
							si++;
						}
					}
					Vector2d mean_sol = sols.rowwise().mean();
					// Convert back to complex
					outputVector[f] = {static_cast<float>(mean_sol[0]), static_cast<float>(mean_sol[1])};
				}
				outputIter.Set(outputVector);
			}
			++inputIter;
			if (m)
				++maskIter;
			++outputIter;
		}
	}

private:
	FirstPassFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
};

class SecondPassFilter : public ImageToImageFilter<VectorImage<complex<float>, 3>, VectorImage<complex<float>, 3>>
{
protected:
	size_t m_flips, m_phases, m_lines = 0;
	bool m_2D = false;

public:
	typedef VectorImage<complex<float>, 3>     TImage;
	typedef Image<float, 3>                    TMask;
	typedef SecondPassFilter                   Self;
	typedef ImageToImageFilter<TImage, TImage> Superclass;
	typedef SmartPointer<Self>                 Pointer;
	typedef typename TImage::RegionType        RegionType;

	itkNewMacro(Self);
	itkTypeMacro(SecondPassFilter, ImageToImageFilter);

	void Set2D(const bool d) { m_2D = d; }
	void SetPhases(const size_t p) {
		if (p < 4)
			throw(runtime_error("Must have a minimum of 4 phase-cycling patterns."));
		if ((p % 2) != 0)
			throw(runtime_error("Number of phases must be even."));
		m_phases = p;
		m_lines = m_phases / 2;
	}
	void SetInput(const TImage *img) { this->SetNthInput(0, const_cast<TImage*>(img)); }
	void SetPass1(const TImage *img) { this->SetNthInput(1, const_cast<TImage*>(img)); }
	void SetMask(const TMask *mask) { this->SetNthInput(2, const_cast<TMask*>(mask)); }
	typename TImage::ConstPointer GetInput() const { return static_cast<const TImage *>(this->ProcessObject::GetInput(0)); }
	typename TImage::ConstPointer GetPass1() const { return static_cast<const TImage *>(this->ProcessObject::GetInput(1)); }
	typename TMask::ConstPointer GetMask() const { return static_cast<const TMask *>(this->ProcessObject::GetInput(2)); }

	virtual void GenerateOutputInformation() override {
		Superclass::GenerateOutputInformation();
		if ((this->GetInput()->GetNumberOfComponentsPerPixel() % m_phases) != 0) {
			throw(std::runtime_error("Input size and number of phases do not match"));
		}
		m_flips = (this->GetInput()->GetNumberOfComponentsPerPixel() / m_phases);
		if (this->GetPass1()->GetNumberOfComponentsPerPixel() != m_flips) {
			throw(std::runtime_error("First passs output has incorrect number of flip-angles"));
		}
		auto op = this->GetOutput();
		op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
		op->SetNumberOfComponentsPerPixel(m_flips);
		op->Allocate();
	}

protected:
	SecondPassFilter() {
		this->SetNumberOfRequiredInputs(2);
		this->SetNumberOfRequiredOutputs(1);
		this->SetNthOutput(0, this->MakeOutput(0));
		this->SetPhases(4);
	}
	~SecondPassFilter() {}

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

	virtual void ThreadedGenerateData(const RegionType &region, ThreadIdType threadId) {
		//std::cout <<  __PRETTY_FUNCTION__ << endl;
		ConstNeighborhoodIterator<TImage>::RadiusType radius;
		radius.Fill(1);
		if (m_2D)
			radius[TImage::ImageDimension - 1] = 0;
		ConstNeighborhoodIterator<TImage> inputIter(radius, this->GetInput(), region);
		ConstNeighborhoodIterator<TImage> pass1Iter(radius, this->GetPass1(), region);

		auto m = this->GetMask();
		ImageRegionConstIterator<TMask> maskIter;
		if (m) {
			maskIter = ImageRegionConstIterator<TMask>(m, region);
		}
		ImageRegionIterator<TImage> outputIter(this->GetOutput(), region);
		while(!inputIter.IsAtEnd()) {
			if (!m || maskIter.Get()) {
				VariableLengthVector<complex<float>> outputVector(m_flips);
				ArrayXd num = ArrayXd::Zero(m_flips), den = ArrayXd::Zero(m_flips);

				for (int p = 0; p < inputIter.Size(); ++p) {
					VariableLengthVector<complex<float>> pass1Vector = pass1Iter.GetPixel(p);
					VariableLengthVector<complex<float>> inputVector = inputIter.GetPixel(p);
					Map<const ArrayXXcf> allInput(inputVector.GetDataPointer(), m_phases, m_flips);
					for (int f = 0; f < m_flips; f++) {
						ArrayXcd a = allInput.col(f).segment(0, m_lines).cast<complex<double>>();
						ArrayXcd b = allInput.col(f).segment(m_lines, m_lines).cast<complex<double>>();
						complex<double> s_i = pass1Vector.GetElement(f);
						for (int li = 0; li < m_lines; li++) {
							complex<double> a_i = a[li];
							complex<double> b_i = b[li];
							num += real(conj(b_i - s_i)*(a_i - b_i) + conj(a_i - b_i)*(b_i - s_i));
							den += real(conj(a_i - b_i)*(a_i - b_i));
						}
					}
				}

				ArrayXd w = -num / (2. * den);
				VariableLengthVector<complex<float>> inputVector = inputIter.GetCenterPixel();
				Map<const ArrayXXcf> allInput(inputVector.GetDataPointer(), m_phases, m_flips);
				for (int f = 0; f < m_flips; f++) {
					if (isfinite(w[f])) {
						ArrayXcd a = allInput.col(f).segment(0, m_lines).cast<complex<double>>();
						ArrayXcd b = allInput.col(f).segment(m_lines, m_lines).cast<complex<double>>();
						for (int li = 0; li < m_lines; li++) {
							complex<double> s = (a[li]*w[f] + (1.f-w[f])*b[li]);
							outputVector[f] += static_cast<complex<float>>(s / static_cast<double>(m_lines));
						}
					}
				}
				outputIter.Set(outputVector);
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
	SecondPassFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
};

} // End namespace itk

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();

	itk::ImageFileReader<itk::Image<float, 3>>::Pointer mask = ITK_NULLPTR;
	auto pass1 = itk::FirstPassFilter::New(); // Needs to be here so save option can be set
	auto pass2 = itk::SecondPassFilter::New(); // Needs to be here so 2D option can be set
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
			case 'F': phaseflip = true; break;
			case 'p': nPhases = atoi(optarg); break;
			case 'a': alternate = true; break;
			case 's':
				switch(*optarg) {
					case 'L': pass1->SetSave(itk::FirstPassFilter::Save::LR); break;
					case 'M': pass1->SetSave(itk::FirstPassFilter::Save::MR); break;
					case 'G': pass1->SetSave(itk::FirstPassFilter::Save::GS); break;
					case 'C': pass1->SetSave(itk::FirstPassFilter::Save::CS); break;
					case 'P': pass1->SetSave(itk::FirstPassFilter::Save::PS); break;
					default:
						cerr << "Unknown desired save image: " << *optarg << endl;
						return EXIT_FAILURE;
				}
				break;
			case '2': do_pass2 = true; if (optarg) pass2->Set2D(true); break;
			case 'T':
				itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg));
				break;
			case 'h':
			case '?': // getopt will print an error message
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
	auto inData = QI::TimeseriesToVectorXF::New();
	auto reorderFlips = QI::ReorderXF::New();
	auto reorderPhase = QI::ReorderXF::New();

	inFile->SetFileName(fname);
	inData->SetInput(inFile->GetOutput());
	reorderFlips->SetInput(inData->GetOutput());       // Does nothing unless stride set
	if (!phaseflip) {
		inData->Update(); // We need to know the vector length to get the number of flips from the number of phases
		size_t nFlips = inData->GetOutput()->GetNumberOfComponentsPerPixel() / nPhases;
		reorderFlips->SetStride(nFlips);
	}
	reorderPhase->SetInput(reorderFlips->GetOutput()); // Does nothing unless stride set
	if (alternate) {
		reorderPhase->SetStride(2);
		reorderPhase->SetBlockSize(nPhases);
	}
	pass1->SetInput(reorderPhase->GetOutput());
	pass1->SetPhases(nPhases);
	pass2->SetPhases(nPhases);
	auto outImage = itk::VectorToImageFilter<complex<float>>::New();
	auto outFile = itk::ImageFileWriter<itk::Image<complex<float>, 4>>::New();
	if (prefix == "")
		prefix = fname.substr(0, fname.find(".nii"));
	string outname = prefix;
	switch (pass1->GetSave()) {
		case itk::FirstPassFilter::Save::LR: outname.append("_lreg"); break;
		case itk::FirstPassFilter::Save::MR: outname.append("_mreg"); break;
		case itk::FirstPassFilter::Save::GS: outname.append("_gs"); break;
		case itk::FirstPassFilter::Save::CS: outname.append("_cs"); break;
		case itk::FirstPassFilter::Save::PS: outname.append("_ps"); break;
	}
	if (do_pass2) {
		outname.append("_2p");
		pass2->SetInput(reorderPhase->GetOutput());
		pass2->SetPass1(pass1->GetOutput());
		outImage->SetInput(pass2->GetOutput());
	} else {
		outImage->SetInput(pass1->GetOutput());
	}
	outname.append(OutExt());
	if (verbose) cout << "Output filename: " << outname << endl;
	outFile->SetInput(outImage->GetOutput());
	outFile->SetFileName(outname);
	if (verbose) cout << "Processing..." << endl;
	outFile->Update();
	if (verbose) cout << "Finished." << endl;
	return EXIT_SUCCESS;
}
