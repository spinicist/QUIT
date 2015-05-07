#ifndef APPLYALGORITHMFILTER_HXX
#define APPLYALGORITHMFILTER_HXX

#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNormalizedCorrelationImageFilter.hxx"

namespace itk {

template<typename TData, typename TAlgo>
ApplyAlgorithmFilter<TData, TAlgo>::ApplyAlgorithmFilter() {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
}

template<typename TData, typename TAlgo>
void ApplyAlgorithmFilter<TData, TAlgo>::Setup() {
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

template<typename TData, typename TAlgo>
void ApplyAlgorithmFilter<TData, TAlgo>::SetDataInput(const size_t i, const TInputImage *image) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_sequence->count()) {
		this->SetNthInput(i, const_cast<TInputImage*>(image));
	} else {
		throw(runtime_error("Data input exceeds range"));
	}
}

template<typename TData, typename TAlgo>
void ApplyAlgorithmFilter<TData, TAlgo>::SetConstInput(const size_t i, const TImage *image) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_algorithm->numConsts()) {
		this->SetNthInput(m_sequence->count() + 1 + i, const_cast<TImage*>(image));
	} else {
		throw(runtime_error("Const input out of range"));
	}
}

template<typename TData, typename TAlgo>
void ApplyAlgorithmFilter<TData, TAlgo>::SetMask(const TImage *image) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	this->SetNthInput(m_sequence->count(), const_cast<TImage*>(image));
}

template<typename TData, typename TAlgo>
auto ApplyAlgorithmFilter<TData, TAlgo>::GetDataInput(const size_t i) const -> typename TInputImage::ConstPointer {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_sequence->count()) {
		return static_cast<const TInputImage *> (this->ProcessObject::GetInput(i));
	} else {
		throw(runtime_error("Get Data Input out of range."));
	}
}

template<typename TData, typename TAlgo>
auto ApplyAlgorithmFilter<TData, TAlgo>::GetConstInput(const size_t i) const -> typename TImage::ConstPointer {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_algorithm->numConsts()) {
		size_t index = m_sequence->count() + 1 + i;
		return static_cast<const TImage *> (this->ProcessObject::GetInput(index));
	} else {
		throw(runtime_error("Get Data Input out of range."));
	}
}

template<typename TData, typename TAlgo>
auto ApplyAlgorithmFilter<TData, TAlgo>::GetMask() const -> typename TImage::ConstPointer {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	return static_cast<const TImage *>(this->ProcessObject::GetInput(m_sequence->count()));
}

template<typename TData, typename TAlgo>
void ApplyAlgorithmFilter<TData, TAlgo>::SetSequence(const shared_ptr<SequenceBase> &seq) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	m_sequence = seq;
}

template<typename TData, typename TAlgo>
void ApplyAlgorithmFilter<TData, TAlgo>::SetAlgorithm(const shared_ptr<TAlgo> &a) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	m_algorithm = a;
}

template<typename TData, typename TAlgo>
DataObject::Pointer ApplyAlgorithmFilter<TData, TAlgo>::MakeOutput(unsigned int idx) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	DataObject::Pointer output;
	if (idx == 0) {
		typename TResidImage::Pointer img = TResidImage::New();
		output = img;
	} else if (idx < (m_algorithm->numOutputs() + 1)) {
		output = (TImage::New()).GetPointer();
	} else {
		std::cerr << "No output " << idx << std::endl;
		output = NULL;
	}
	return output.GetPointer();
}

template<typename TData, typename TAlgo>
auto ApplyAlgorithmFilter<TData, TAlgo>::GetOutput(const size_t i) -> TImage *{
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_algorithm->numOutputs()) {
		return dynamic_cast<TImage *>(this->ProcessObject::GetOutput(i+1));
	} else {
		throw(runtime_error("Requested output " + to_string(i) + " is past maximum (" + to_string(m_algorithm->numOutputs()) + ")"));
	}
}

template<typename TData, typename TAlgo>
auto ApplyAlgorithmFilter<TData, TAlgo>::GetResidOutput() -> TResidImage *{
	return dynamic_cast<TResidImage *>(this->ProcessObject::GetOutput(0));
}

template<typename TData, typename TAlgo>
void ApplyAlgorithmFilter<TData, TAlgo>::Update() {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	Superclass::Update();
}

template<typename TData, typename TAlgo>
void ApplyAlgorithmFilter<TData, TAlgo>::GenerateOutputInformation() {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	Superclass::GenerateOutputInformation();
	size_t size = 0;
	for (size_t i = 0; i < m_sequence->count(); i++) {
		size += this->GetDataInput(i)->GetNumberOfComponentsPerPixel();
	}
	if (m_sequence->size() != size) {
		throw(std::runtime_error("Sequence size (" + to_string(m_sequence->size()) + ") does not match input size (" + to_string(size) + ")"));
	}

	for (size_t i = 0; i < (m_algorithm->numOutputs()); i++) {
		auto op = this->GetOutput(i);
		op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
		op->Allocate();
	}
	auto r = this->GetResidOutput();
	r->SetRegions(this->GetInput()->GetLargestPossibleRegion());
	r->SetNumberOfComponentsPerPixel(size);
	r->Allocate();

}

template<typename TData, typename TAlgo>
void ApplyAlgorithmFilter<TData, TAlgo>::ThreadedGenerateData(const RegionType & region, ThreadIdType threadId) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	vector<ImageRegionConstIterator<TInputImage>> dataIters(m_sequence->count());
	for (size_t i = 0; i < m_sequence->count(); i++) {
		dataIters[i] = ImageRegionConstIterator<TInputImage>(this->GetDataInput(i), region);
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
	ImageRegionIterator<TResidImage> residIter(this->GetResidOutput(), region);

	VectorXd constants = m_algorithm->defaultConsts();
	while(!dataIters[0].IsAtEnd()) {
		typename TImage::ConstPointer m = this->GetMask();
		if (!m || maskIter.Get()) {
			for (size_t i = 0; i < constIters.size(); i++) {
				if (this->GetConstInput(i))
					constants[i] = constIters[i].Get();
			}

			typename TAlgo::TInputVector allData(m_sequence->size());
			typedef typename TAlgo::TInput TInput;
			size_t dataIndex = 0;
			for (size_t i = 0; i < m_sequence->count(); i++) {
				VariableLengthVector<TData> dataVector = dataIters[i].Get();
				Map<const Eigen::Array<TData, Eigen::Dynamic, 1>> data(dataVector.GetDataPointer(), dataVector.Size());
				allData.segment(dataIndex, data.rows()) = data.template cast<TInput>();
			}
			VectorXd outputs(m_algorithm->numOutputs());

			ArrayXd resids(m_sequence->size());

			m_algorithm->apply(m_sequence, allData, constants, outputs, resids);

			for (size_t i = 0; i < m_algorithm->numOutputs(); i++) {
				outputIters[i].Set(static_cast<float>(outputs[i]));
			}
			ArrayXf residF = resids.cast<float>();
			VariableLengthVector<float> residVector(residF.data(), m_sequence->size());

			residIter.Set(residVector);
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
		++residIter;
	}
}
} // namespace ITK

#endif // APPLYALGORITHMFILTER_HXX

