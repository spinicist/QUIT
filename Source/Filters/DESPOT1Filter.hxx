#ifndef DESPOT1FILTER_HXX
#define DESPOT1FILTER_HXX

#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNormalizedCorrelationImageFilter.hxx"

namespace itk {

template<typename TVectorImage, typename TImage>
DESPOT1Filter<TVectorImage, TImage>::DESPOT1Filter()
{
	//std::cout << __PRETTY_FUNCTION__ << endl;
	this->SetNumberOfRequiredInputs(1);
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
void DESPOT1Filter<TVectorImage, TImage>::SetSequence(const shared_ptr<SPGRSimple> &seq) {
	//std::cout << __PRETTY_FUNCTION__ << endl;
	m_sequence = seq;
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetAlgorithm(const shared_ptr<Algorithm> &a) {
	//std::cout << __PRETTY_FUNCTION__ << endl;
	m_algorithm = a;
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
	if (m_sequence->size() != size) {
		throw(std::runtime_error("Specified number of flip-angles does not match number of volumes in input."));
	}

	for (size_t i = 0; i < 3; i++) {
		this->GetOutput(i)->SetRegions(this->GetInput()->GetLargestPossibleRegion());
		this->GetOutput(i)->Allocate();
	}
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::ThreadedGenerateData(const RegionType & region, ThreadIdType threadId) {
	typename TVectorImage::ConstPointer spgrData = this->GetInput();

	typename TImage::ConstPointer maskData = this->GetMask();
	typename TImage::ConstPointer B1Data   = this->GetB1();

	typename TImage::Pointer T1Data = this->GetOutput(0);
	typename TImage::Pointer PDData = this->GetOutput(1);
	typename TImage::Pointer ResData = this->GetOutput(2);

	ImageRegionConstIterator<TVectorImage> SPGRIter(spgrData, region);
	ImageRegionConstIterator<TImage> maskIter, B1Iter;
	if (maskData)
		maskIter = ImageRegionConstIterator<TImage>(maskData, region);
	if (B1Data)
		B1Iter = ImageRegionConstIterator<TImage>(B1Data, region);
	ImageRegionIterator<TImage> T1Iter(T1Data, region);
	ImageRegionIterator<TImage> PDIter(PDData, region);
	ImageRegionIterator<TImage> ResIter(ResData, region);

	while(!T1Iter.IsAtEnd()) {
		if (!maskData || maskIter.Get()) {
			VectorXd consts = VectorXd::Ones(1);
			if (B1Data)
				consts[0] = B1Iter.Get();

			VariableLengthVector<float> signalVector = SPGRIter.Get();
			Map<const ArrayXf> signalf(signalVector.GetDataPointer(), m_sequence->size());

			VectorXd outputs(m_algorithm->numOutputs());
			ArrayXd resids(m_sequence->size());

			m_algorithm->apply(m_sequence, signalf.cast<double>(), consts, outputs, resids);

			/*if (all_residuals) {
				ResidsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = resids.cast<float>();
			}*/
			T1Iter.Set(static_cast<float>(outputs[1]));
			PDIter.Set(static_cast<float>(outputs[0]));
			ResIter.Set(static_cast<float>(sqrt(resids.square().sum() / resids.rows())));
		} else {
			//T1Iter.Set(0);
			//PDIter.Set(0);
			//ResIter.Set(0);
		}
		++SPGRIter;
		if (maskData)
			++maskIter;
		if (B1Data)
			++B1Iter;
		++T1Iter; ++PDIter; ++ResIter;
	}
}
} // namespace ITK

#endif // DESPOT1FILTER_HXX

