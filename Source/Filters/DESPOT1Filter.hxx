#ifndef DESPOT1FILTER_HXX
#define DESPOT1FILTER_HXX

#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNormalizedCorrelationImageFilter.hxx"

namespace itk {

template<typename TVectorImage, typename TImage>
DESPOT1Filter<TVectorImage, TImage>::DESPOT1Filter() {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::Setup() {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (!(m_sequence && m_algorithm))
		throw(runtime_error("Sequence and Algorithm must be set first"));

	// +1 is for mask.
	// Inputs go: Data 0, Data 1, ..., Mask, Const 0, Const 1, ...
	size_t totalInputs = m_sequence->count() + m_algorithm->numConsts() + 1;
	this->SetNumberOfRequiredInputs(m_sequence->count());
	// +1 is for residuals vector
	size_t totalOutputs = m_algorithm->numOutputs() + 1;
	this->SetNumberOfRequiredOutputs(totalOutputs);

	for (size_t i = 0; i < totalOutputs; i++) {
		this->SetNthOutput(i, this->MakeOutput(i));
	}
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetDataInput(const size_t i, const TVectorImage *image) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_sequence->count()) {
		this->SetNthInput(i, const_cast<TVectorImage*>(image));
	} else {
		throw(runtime_error("Data input exceeds range"));
	}
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetConstInput(const size_t i, const TImage *image) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_algorithm->numConsts()) {
		this->SetNthInput(m_sequence->count() + 1 + i, const_cast<TImage*>(image));
	} else {
		throw(runtime_error("Const input out of range"));
	}
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetMask(const TImage *image) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	this->SetNthInput(m_sequence->count(), const_cast<TImage*>(image));
}

template< typename TVectorImage, typename TImage>
typename TVectorImage::ConstPointer DESPOT1Filter<TVectorImage, TImage>::GetDataInput(const size_t i) const {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_sequence->count()) {
		return static_cast<const TVectorImage *> (this->ProcessObject::GetInput(i));
	} else {
		throw(runtime_error("Get Data Input out of range."));
	}
}

template< typename TVectorImage, typename TImage>
typename TImage::ConstPointer DESPOT1Filter<TVectorImage, TImage>::GetConstInput(const size_t i) const {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_algorithm->numConsts()) {
		size_t index = m_sequence->count() + 1 + i;
		return static_cast<const TImage *> (this->ProcessObject::GetInput(index));
	} else {
		throw(runtime_error("Get Data Input out of range."));
	}
}

template< typename TVectorImage, typename TImage>
typename TImage::ConstPointer DESPOT1Filter<TVectorImage, TImage>::GetMask() const {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	return static_cast<const TImage *>(this->ProcessObject::GetInput(m_sequence->count()));
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetSequence(const shared_ptr<SPGRSimple> &seq) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	m_sequence = seq;
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::SetAlgorithm(const shared_ptr<Algorithm> &a) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	m_algorithm = a;
}

template<typename TVectorImage, typename TImage>
DataObject::Pointer DESPOT1Filter<TVectorImage, TImage>::MakeOutput(unsigned int idx) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	cout << "idx = " << idx << endl;
	DataObject::Pointer output;
	if (idx == 0) {
		typename TVectorImage::Pointer img = TVectorImage::New();
		img->SetNumberOfComponentsPerPixel(m_sequence->size());
		output = img;
	} else if (idx < (m_algorithm->numOutputs() + 1)) {
		output = (TImage::New()).GetPointer();
	} else {
		std::cerr << "No output " << idx << std::endl;
		output = NULL;
	}
	return output.GetPointer();
}

template< typename TVectorImage, typename TImage>
TImage *DESPOT1Filter<TVectorImage, TImage>::GetOutput(const size_t i) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_algorithm->numOutputs()) {
		return dynamic_cast<TImage *>(this->ProcessObject::GetOutput(i+1) );
	} else {
		throw(runtime_error("Requested output " + to_string(i) + " is past maximum (" + to_string(m_algorithm->numOutputs()) + ")"));
	}
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::Update() {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	//std::cout <<  static_cast<const TVectorImage *>(this->ProcessObject::GetInput(0))->GetDirection() << endl;
	//std::cout <<  static_cast<const TImage *>(this->ProcessObject::GetInput(1))->GetDirection() << endl;
	//std::cout <<  static_cast<const TImage *>(this->ProcessObject::GetInput(2))->GetDirection() << endl;
	Superclass::Update();
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::GenerateOutputInformation() {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	Superclass::GenerateOutputInformation();
	size_t size = 0;
	for (size_t i = 0; i < m_sequence->count(); i++) {
		size += this->GetDataInput(i)->GetNumberOfComponentsPerPixel();
	}
	if (m_sequence->size() != size) {
		throw(std::runtime_error("Specified number of flip-angles does not match number of volumes in input."));
	}

	for (size_t i = 0; i < (m_algorithm->numOutputs()); i++) {
		const auto op = this->GetOutput(i);
		op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
		op->Allocate();
	}
}

template<typename TVectorImage, typename TImage>
void DESPOT1Filter<TVectorImage, TImage>::ThreadedGenerateData(const RegionType & region, ThreadIdType threadId) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	vector<ImageRegionConstIterator<TVectorImage>> dataIters(m_sequence->count());
	for (size_t i = 0; i < m_sequence->count(); i++) {
		dataIters[i] = ImageRegionConstIterator<TVectorImage>(this->GetDataInput(i), region);
	}
	ImageRegionConstIterator<TImage> maskIter;
	if (this->GetMask()) {
		maskIter = ImageRegionConstIterator<TImage>(this->GetMask(), region);
	}
	vector<ImageRegionConstIterator<TImage>> constIters(m_algorithm->numConsts());
	for (size_t i = 0; i < m_algorithm->numConsts(); i++) {
		typename TImage::ConstPointer c = this->GetConstInput(i);
		if (c) {
			constIters[i] = ImageRegionConstIterator<TImage>(c, region);
		}
	}
	vector<ImageRegionIterator<TImage>> outputIters(m_algorithm->numOutputs());
	for (size_t i = 0; i < m_algorithm->numOutputs(); i++) {
		outputIters[i] = ImageRegionIterator<TImage>(this->GetOutput(i), region);
	}

	VectorXd constants = m_algorithm->defaultConsts();
	while(!dataIters[0].IsAtEnd()) {
		typename TImage::ConstPointer m = this->GetMask();
		if (!m || maskIter.Get()) {
			for (size_t i = 0; i < constIters.size(); i++) {
				if (this->GetConstInput(i))
					constants[i] = constIters[i].Get();
			}

			ArrayXf allData(m_sequence->size());
			size_t dataIndex = 0;
			for (size_t i = 0; i < m_sequence->count(); i++) {
				VariableLengthVector<float> dataVector = dataIters[i].Get();
				Map<const ArrayXf> data(dataVector.GetDataPointer(), dataVector.Size());
				allData.segment(dataIndex, data.rows()) = data;
			}
			VectorXd outputs(m_algorithm->numOutputs());
			ArrayXd resids(m_sequence->size());

			m_algorithm->apply(m_sequence, allData.cast<double>(), constants, outputs, resids);

			/*if (all_residuals) {
				ResidsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = resids.cast<float>();
			}*/

			for (size_t i = 0; i < m_algorithm->numOutputs(); i++) {
				outputIters[i].Set(static_cast<float>(outputs[i]));
				//ResIter.Set(static_cast<float>(sqrt(resids.square().sum() / resids.rows())));
			}
		} else {
			//T1Iter.Set(0);
			//PDIter.Set(0);
			//ResIter.Set(0);
		}
		for (size_t i = 0; i < m_sequence->count(); i++) {
			++dataIters[i];
		}
		if (this->GetMask())
			++maskIter;
		for (size_t i = 0; i < m_algorithm->numConsts(); i++) {
			if (this->GetConstInput(i))
				++constIters[i];
		}
		for (size_t i = 0; i < m_algorithm->numOutputs(); i++) {
			++outputIters[i];
		}
	}
}
} // namespace ITK

#endif // DESPOT1FILTER_HXX

