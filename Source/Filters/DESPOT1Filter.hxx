#ifndef DESPOT1FILTER_HXX
#define DESPOT1FILTER_HXX

#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk {

template<typename TVectorImage, typename TImage>
DESPOT1Filter<TVectorImage, TImage>::DESPOT1Filter() :
	m_sequence({},0)
{
	//std::cout << __PRETTY_FUNCTION__ << endl;
	this->SetNumberOfRequiredInputs(3);
	this->SetNumberOfRequiredOutputs(3);

	this->SetNthOutput(0, this->MakeOutput(0));
	this->SetNthOutput(1, this->MakeOutput(1));
	this->SetNthOutput(2, this->MakeOutput(2));
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetInput(const TVectorImage *image) {
	//std::cout << __PRETTY_FUNCTION__ << endl;
	this->SetNthInput(0, const_cast<TVectorImage*>(image));
	//std::cout << this->GetInput()->GetDirection() << endl;
}
template< typename TVectorImage, typename TImage>
typename TVectorImage::ConstPointer DESPOT1Filter<TVectorImage, TImage>::GetInput() {
	//std::cout << __PRETTY_FUNCTION__ << endl;
	return static_cast<const TVectorImage *> (this->ProcessObject::GetInput(0));
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetMask(const TImage* image) {
	//std::cout << __PRETTY_FUNCTION__ << endl;
	this->SetNthInput(1, const_cast<TImage*>(image));
	//std::cout << this->GetMask()->GetDirection() << endl;
}
template< typename TVectorImage, typename TImage>
typename TImage::ConstPointer DESPOT1Filter<TVectorImage, TImage>::GetMask() {
	//std::cout << __PRETTY_FUNCTION__ << endl;
	return static_cast<const TImage *>(this->ProcessObject::GetInput(1));
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetB1(const TImage* image) {
	//std::cout << __PRETTY_FUNCTION__ << endl;
	this->SetNthInput(2, const_cast<TImage*>(image));
	//std::cout << this->GetB1()->GetDirection() << endl;
}
template<typename TVectorImage, typename TImage>
typename TImage::ConstPointer DESPOT1Filter<TVectorImage, TImage>::GetB1() {
	//std::cout << __PRETTY_FUNCTION__ << endl;
	return static_cast<const TImage *>(this->ProcessObject::GetInput(2));
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetSequence(const SPGRSimple &seq) {
	//std::cout << __PRETTY_FUNCTION__ << endl;
	m_sequence = seq;
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetAlgorithm(const Algos &a) {
	//std::cout << __PRETTY_FUNCTION__ << endl;
	m_algorithm = a;
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetIterations(const size_t &n) {
	//std::cout << __PRETTY_FUNCTION__ << endl;
	m_iterations = n;
}

template<typename TVectorImage, typename TImage>
DataObject::Pointer DESPOT1Filter<TVectorImage, TImage>::MakeOutput(unsigned int idx) {
	//std::cout << __PRETTY_FUNCTION__ << endl;
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
TImage* DESPOT1Filter<TVectorImage, TImage>::GetOutput(const size_t i) {
	//std::cout << __PRETTY_FUNCTION__ << endl;
	return dynamic_cast<TImage *>(this->ProcessObject::GetOutput(i) );
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::Update() {
	std::cout << __PRETTY_FUNCTION__ << endl;
	//std::cout << static_cast<const TVectorImage *>(this->ProcessObject::GetInput(0))->GetDirection() << endl;
	//std::cout << static_cast<const TImage *>(this->ProcessObject::GetInput(1))->GetDirection() << endl;
	//std::cout << static_cast<const TImage *>(this->ProcessObject::GetInput(2))->GetDirection() << endl;
	Superclass::Update();
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::GenerateOutputInformation() {
	Superclass::GenerateOutputInformation();
	auto size = this->GetInput()->GetNumberOfComponentsPerPixel();
	cout << "Size check: " << m_sequence.size() << "/" << size << endl;
	if (m_sequence.size() != size) {
		throw(std::runtime_error("Specified number of flip-angles does not match number of volumes in input."));
	}
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::GenerateData() {
	//std::cout << __PRETTY_FUNCTION__ << endl;
	typename TVectorImage::ConstPointer spgrData = this->GetInput();

	typename TImage::ConstPointer maskData = this->GetMask();
	typename TImage::ConstPointer B1Data   = this->GetB1();

	typename TImage::Pointer T1Data = this->GetOutput(0);
	typename TImage::Pointer PDData = this->GetOutput(1);
	typename TImage::Pointer ResData = this->GetOutput(2);

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

#endif // DESPOT1FILTER_HXX

