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

#include "itkImageFileReader.h"

#include "Util.h"
#include "Filters/ImageToVectorFilter.h"
#include "Filters/VectorToImageFilter.h"

using namespace std;
using namespace Eigen;
using namespace QUITK;

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

static bool verbose = false, pass2 = false, pass22d = false, alternate = false;
static string prefix;
const struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"flip", required_argument, 0, 'F'},
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
	bool m_phase_first, m_alternate = false;
	Save m_mode = Save::LR;
public:
	/** Standard class typedefs. */
	typedef VectorImage<complex<float>, 3>      TImage;
	typedef Image<float, 3>                     TMask;
	typedef FirstPassFilter                     Self;
	typedef ImageToImageFilter<TImage, TImage> Superclass;
	typedef SmartPointer<Self>                  Pointer;
	typedef typename TImage::RegionType         RegionType;

	itkNewMacro(Self); /** Method for creation through the object factory. */
	itkTypeMacro(FirstPassFilter, ImageToImageFilter); /** Run-time type information (and related methods). */

	void SetSave(Save s) { m_mode = s; }
	Save GetSave() { return m_mode; }
	void SetPhaseFirst(const bool f) { m_phase_first = f; }
	void SetAlternate(const bool a) { m_alternate = a; }
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

				ArrayXXcf allInput;
				if (m_phase_first) {
					allInput = Eigen::Map<const Eigen::Array<complex<float>, Eigen::Dynamic, Eigen::Dynamic>>(inputVector.GetDataPointer(), m_phases, m_flips);
				} else {
					allInput = Eigen::Map<const Eigen::Array<complex<float>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(inputVector.GetDataPointer(), m_phases, m_flips);
				}
				for (int f = m_flips - 1; f > -1; f--) {
					ArrayXcf thisFlip = allInput.col(f);
					ArrayXcd a, b;
					if (m_alternate) {
						a = Map<Eigen::ArrayXcf, Eigen::Unaligned, Eigen::Stride<2, 1>>(thisFlip.data(), m_lines).cast<complex<double>>();
						b = Map<Eigen::ArrayXcf, Eigen::Unaligned, Eigen::Stride<2, 1>>(thisFlip.data()+1, m_lines).cast<complex<double>>();
					} else {
						a = thisFlip.segment(0, m_lines).cast<complex<double>>();
						b = thisFlip.segment(m_lines, m_lines).cast<complex<double>>();
					}
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

} // End namespace itk

/*
class Pass2 : public Algorithm {
		if (pass2) {
			for (size_t vk = 0; vk < d[2]; vk++) {
			//size_t vk = 27;
				function<void (const size_t, const size_t)> processVox = [&] (const size_t vi, const size_t vj) {
					complex<float> sp(0.,0.);
					if (!maskFile || (maskData[{vi,vj,vk}])) {
						for (size_t li = 0; li < nLines; li++) {
							float num = 0, den = 0;
							for (int k = ((!pass22d && (vk > 0)) ? -1 : 0); k < ((!pass22d && (vk < d[2] - 1)) ? 2 : 1); k++) {
								for (int j = (vj > 0 ? -1 : 0); j < (vj < d[1] - 1 ? 2 : 1); j++) {
									for (int i = (vi > 0 ? -1 : 0); i < (vi < d[0] - 1 ? 2 : 1); i++) {
										idx_t idx; idx << vi + i,vj + j,vk + k,0,0;
										idx[flip_dim] = vol;
										idx[phase_dim] = li;
										complex<float> a_i = aData[idx];
										complex<float> b_i = bData[idx];
										complex<float> s_i = output[{vi+i,vj+j,vk+k,vol}];
										num += real(conj(b_i - s_i)*(a_i - b_i) + conj(a_i - b_i)*(b_i - s_i));
										den += real(conj(a_i - b_i)*(a_i - b_i));
									}
								}
							}
							float w = -num / (2. * den);
							if (isfinite(w)) {
								idx_t sidx; sidx << vi,vj,vk,0,0;
								sidx[flip_dim] = vol;
								sidx[phase_dim] = li;
								complex<float> s = (aData[sidx]*w + (1.f-w)*bData[sidx]);
								sp += (s / static_cast<float>(nLines));
							}
						}
					}
					second_pass[{vi,vj,vk,vol}] = sp;
				};
				threads.for_loop2(processVox, 0, d[0], 1, 0, d[1], 1);
				if (threads.interrupted())
					break;
			}
		}

};*/


//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();

	itk::ImageFileReader<itk::Image<float, 3>>::Pointer mask = ITK_NULLPTR;
	auto pass1 = itk::FirstPassFilter::New();
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'm':
				if (verbose) cout << "Reading mask file " << optarg << endl;
				mask = itk::ImageFileReader<itk::Image<float, 3>>::New();
				mask->SetFileName(optarg);
				break;
			case 'o':
				prefix = optarg;
				cout << "Output prefix will be: " << prefix << endl;
				break;
			case 'F':
				pass1->SetPhaseFirst(true); break;
			case 'p':
				pass1->SetPhases(atoi(optarg));
				break;
			case 'a': pass1->SetAlternate(true); break;
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
			case '2': pass2 = true; if (optarg) pass22d = true; break;
			case 'T':
				itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg));
				break;
			case 'h':
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
		}
	}
	if (pass22d) cout << "2D second pass selected." << endl;
	if ((argc - optind) != 1) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		return EXIT_FAILURE;
	}
	if (verbose) cout << "Opening input file: " << argv[optind] << endl;
	string fname(argv[optind++]);

	auto inFile = itk::ImageFileReader<itk::Image<complex<float>, 4>>::New();
	auto inImage = ImageToVectorFilter<complex<float>>::New();
	inFile->SetFileName(fname);
	inImage->SetInput(inFile->GetOutput());
	pass1->SetInput(inImage->GetOutput());
	auto outImage = VectorToImageFilter<complex<float>>::New();
	auto outFile = itk::ImageFileWriter<itk::Image<complex<float>, 4>>::New();
	outImage->SetInput(pass1->GetOutput());
	outFile->SetInput(outImage->GetOutput());
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
	if (pass2) outname.append("_2p");
	outname.append(OutExt());
	if (verbose) cout << "Output filename: " << outname << endl;
	outFile->SetFileName(outname);
	outFile->Update();
	if (verbose) cout << "Finished." << endl;
	return EXIT_SUCCESS;
}
