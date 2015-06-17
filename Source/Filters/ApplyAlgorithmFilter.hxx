#ifndef APPLYALGORITHMFILTER_HXX
#define APPLYALGORITHMFILTER_HXX

#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkProgressReporter.h"

namespace itk {

template<typename TVImage, typename TAlgo>
ApplyAlgorithmFilter<TVImage, TAlgo>::ApplyAlgorithmFilter() {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
}

template<typename TVImage, typename TAlgo>
void ApplyAlgorithmFilter<TVImage, TAlgo>::SetAlgorithm(const shared_ptr<TAlgo> &a) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	m_algorithm = a;
	// +1 is for mask.
	// Inputs go: Data 0, Data 1, ..., Mask, Const 0, Const 1, ...
	this->SetNumberOfRequiredInputs(a->numInputs());
	// +1 is for residuals vector
	size_t totalOutputs = m_algorithm->numOutputs() + 1;
	this->SetNumberOfRequiredOutputs(totalOutputs);
	for (size_t i = 0; i < totalOutputs; i++) {
		this->SetNthOutput(i, this->MakeOutput(i));
	}
}

template<typename TVImage, typename TAlgo>
shared_ptr<const TAlgo> ApplyAlgorithmFilter<TVImage, TAlgo>::GetAlgorithm() const { return m_algorithm; }

template<typename TVImage, typename TAlgo>
void ApplyAlgorithmFilter<TVImage, TAlgo>::SetScaleToMean(const bool s) { m_scale_to_mean = s; }

template<typename TVImage, typename TAlgo>
bool ApplyAlgorithmFilter<TVImage, TAlgo>::GetScaleToMean() const { return m_scale_to_mean; }

template<typename TVImage, typename TAlgo>
void ApplyAlgorithmFilter<TVImage, TAlgo>::SetInput(const size_t i, const TVImage *image) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_algorithm->numInputs()) {
		this->SetNthInput(i, const_cast<TVectorImage*>(image));
	} else {
		throw(runtime_error("Data input exceeds range"));
	}
}

template<typename TVImage, typename TAlgo>
void ApplyAlgorithmFilter<TVImage, TAlgo>::SetConst(const size_t i, const TImage *image) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_algorithm->numConsts()) {
		this->SetNthInput(m_algorithm->numInputs() + 1 + i, const_cast<TImage*>(image));
	} else {
		throw(runtime_error("ConstInput " + to_string(i) + " out of range (there are " + to_string(m_algorithm->numConsts()) + " inputs)"));
	}
}

template<typename TVImage, typename TAlgo>
void ApplyAlgorithmFilter<TVImage, TAlgo>::SetMask(const TImage *image) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	this->SetNthInput(m_algorithm->numInputs(), const_cast<TImage*>(image));
}

template<typename TVImage, typename TAlgo>
void ApplyAlgorithmFilter<TVImage, TAlgo>::SetSlices(const int start, const int stop) {
	m_startSlice = start;
	m_stopSlice = stop;
}

template<typename TVImage, typename TAlgo>
auto ApplyAlgorithmFilter<TVImage, TAlgo>::GetInput(const size_t i) const -> typename TVImage::ConstPointer {
	////std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_algorithm->numInputs()) {
		return static_cast<const TVImage *> (this->ProcessObject::GetInput(i));
	} else {
		throw(runtime_error(__PRETTY_FUNCTION__ +
		                    string("Input ") + to_string(i) +
		                    " out of range (" + to_string(m_algorithm->numInputs()) + ")"));
	}
}

template<typename TVImage, typename TAlgo>
auto ApplyAlgorithmFilter<TVImage, TAlgo>::GetConst(const size_t i) const -> typename TImage::ConstPointer {
	////std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_algorithm->numConsts()) {
		size_t index = m_algorithm->numInputs() + 1 + i;
		return static_cast<const TImage *> (this->ProcessObject::GetInput(index));
	} else {
		throw(runtime_error("Get Data Input out of range."));
	}
}

template<typename TVImage, typename TAlgo>
auto ApplyAlgorithmFilter<TVImage, TAlgo>::GetMask() const -> typename TImage::ConstPointer {
	////std::cout <<  __PRETTY_FUNCTION__ << endl;
	return static_cast<const TImage *>(this->ProcessObject::GetInput(m_algorithm->numInputs()));
}

template<typename TVImage, typename TAlgo>
DataObject::Pointer ApplyAlgorithmFilter<TVImage, TAlgo>::MakeOutput(unsigned int idx) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	DataObject::Pointer output;
	if (idx == 0) {
		auto img = TVImage::New();
		output = img;
	} else if (idx < (m_algorithm->numOutputs() + 1)) {
		output = (TImage::New()).GetPointer();
	} else {
		std::cerr << "No output " << idx << std::endl;
		output = NULL;
	}
	return output.GetPointer();
}

template<typename TVImage, typename TAlgo>
auto ApplyAlgorithmFilter<TVImage, TAlgo>::GetOutput(const size_t i) -> TImage *{
	////std::cout <<  __PRETTY_FUNCTION__ << endl;
	if (i < m_algorithm->numOutputs()) {
		return dynamic_cast<TImage *>(this->ProcessObject::GetOutput(i+1));
	} else {
		throw(runtime_error("Requested output " + to_string(i) + " is past maximum (" + to_string(m_algorithm->numOutputs()) + ")"));
	}
}

template<typename TVImage, typename TAlgo>
auto ApplyAlgorithmFilter<TVImage, TAlgo>::GetResidOutput() -> TVImage *{
	return dynamic_cast<TVImage *>(this->ProcessObject::GetOutput(0));
}

template<typename TVImage, typename TAlgo>
void ApplyAlgorithmFilter<TVImage, TAlgo>::GenerateOutputInformation() {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	Superclass::GenerateOutputInformation();
	size_t size = 0;
	for (size_t i = 0; i < m_algorithm->numInputs(); i++) {
		size += this->GetInput(i)->GetNumberOfComponentsPerPixel();
	}
	if (m_algorithm->dataSize() != size) {
		throw(std::runtime_error("Sequence size (" + to_string(m_algorithm->dataSize()) + ") does not match input size (" + to_string(size) + ")"));
	}

	auto region = this->GetInput(0)->GetLargestPossibleRegion();
	const int LastDim = ImageDimension - 1;
	region.GetModifiableIndex()[LastDim] = m_startSlice;
	if (m_stopSlice != 0)
		region.GetModifiableSize()[LastDim] = m_stopSlice - m_startSlice;
	else
		region.GetModifiableSize()[LastDim] = region.GetSize()[LastDim] - m_startSlice;
	for (size_t i = 0; i < m_algorithm->numOutputs(); i++) {
		auto op = this->GetOutput(i);
		op->SetRegions(region);
		//std::cout << "ALLOCATING OUTPUT " << i << std::endl;
		//std::cout << region << std::endl;
		op->Allocate();
	}
	auto r = this->GetResidOutput();
	r->SetRegions(region);
	r->SetNumberOfComponentsPerPixel(size);
	//std::cout << "ALLOCATING RESID" << std::endl;
	//std::cout << region << std::endl;
	//std::cout << size << std::endl;
	r->Allocate();
	//std::cout <<  "Finished " << __PRETTY_FUNCTION__ << endl;
}

template<typename TVImage, typename TAlgo>
void ApplyAlgorithmFilter<TVImage, TAlgo>::ThreadedGenerateData(const TRegion &region, ThreadIdType threadId) {
	//std::cout <<  __PRETTY_FUNCTION__ << std::endl;
	//std::cout << "Thread " << threadId << std::endl;
	//std::cout << region << std::endl;

	ProgressReporter progress(this, threadId, region.GetNumberOfPixels(), 10);

	vector<ImageRegionConstIterator<TVImage>> dataIters(m_algorithm->numInputs());
	for (size_t i = 0; i < m_algorithm->numInputs(); i++) {
		dataIters[i] = ImageRegionConstIterator<TVImage>(this->GetInput(i), region);
	}

	ImageRegionConstIterator<TImage> maskIter;
	const auto mask = this->GetMask();
	if (mask) {
		maskIter = ImageRegionConstIterator<TImage>(mask, region);
	}
	vector<ImageRegionConstIterator<TImage>> constIters(m_algorithm->numConsts());
	for (size_t i = 0; i < m_algorithm->numConsts(); i++) {
		typename TImage::ConstPointer c = this->GetConst(i);
		if (c) {
			constIters[i] = ImageRegionConstIterator<TImage>(c, region);
		}
	}
	vector<ImageRegionIterator<TImage>> outputIters(m_algorithm->numOutputs());
	for (size_t i = 0; i < m_algorithm->numOutputs(); i++) {
		outputIters[i] = ImageRegionIterator<TImage>(this->GetOutput(i), region);
	}
	ImageRegionIterator<TVImage> residIter(this->GetResidOutput(), region);

	typedef typename TAlgo::TArray TArray;
	while(!dataIters[0].IsAtEnd()) {
		TArray outputs = TArray::Zero(m_algorithm->numOutputs());
		TArray resids =  TArray::Zero(m_algorithm->dataSize());
		if (!mask || maskIter.Get()) {
			TArray constants = m_algorithm->defaultConsts();
			for (size_t i = 0; i < constIters.size(); i++) {
				if (this->GetConst(i)) {
					constants[i] = constIters[i].Get();
				}
			}
			typename TAlgo::TInput allData(m_algorithm->dataSize());
			size_t dataIndex = 0;
			for (size_t i = 0; i < m_algorithm->numInputs(); i++) {
				VariableLengthVector<TPixel> dataVector = dataIters[i].Get();
				Map<const Eigen::Array<TPixel, Eigen::Dynamic, 1>> data(dataVector.GetDataPointer(), dataVector.Size());
				Eigen::Array<TPixel, Eigen::Dynamic, 1> scaled = data;
				if (m_scale_to_mean) {
					scaled /= scaled.mean();
				}
				allData.segment(dataIndex, data.rows()) = scaled.template cast<typename TAlgo::TScalar>();
				dataIndex += data.rows();
			}
			m_algorithm->apply(allData, constants, outputs, resids);
		}
		if (this->GetMask())
			++maskIter;
		for (size_t i = 0; i < m_algorithm->numInputs(); i++) {
			++dataIters[i];
		}
		for (size_t i = 0; i < m_algorithm->numConsts(); i++) {
			if (this->GetConst(i))
				++constIters[i];
		}
		for (size_t i = 0; i < m_algorithm->numOutputs(); i++) {
			outputIters[i].Set(static_cast<float>(outputs[i]));
			++outputIters[i];
		}
		ArrayXf residF = resids.template cast<float>();
		VariableLengthVector<float> residVector(residF.data(), m_algorithm->dataSize());
		residIter.Set(residVector);
		++residIter;
		progress.CompletedPixel();
	}
	//std::cout << "Finished " << std::endl;
}
} // namespace ITK

#endif // APPLYALGORITHMFILTER_HXX
