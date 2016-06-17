#ifndef APPLYALGORITHMFILTER_HXX
#define APPLYALGORITHMFILTER_HXX

#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkProgressReporter.h"

namespace itk {

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::ApplyAlgorithmFilter() {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::SetAlgorithm(const std::shared_ptr<TAlgorithm> &a) {
	//std::cout <<  __PRETTY_FUNCTION__ << endl;
	m_algorithm = a;
	// Inputs go: Data 0, Data 1, ..., Mask, Const 0, Const 1, ...
    // Only the data inputs are required, the others are optional
    this->SetNumberOfRequiredInputs(a->numInputs());
    // Outputs go: AllResiduals, Residual, Iterations, Parameter 0, Parameter 1, ...
    size_t totalOutputs = StartOutputs + m_algorithm->numOutputs();
	this->SetNumberOfRequiredOutputs(totalOutputs);
	for (size_t i = 0; i < totalOutputs; i++) {
		this->SetNthOutput(i, this->MakeOutput(i));
	}
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
std::shared_ptr<const TAlgorithm> ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::GetAlgorithm() const { return m_algorithm; }

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::SetScaleToMean(const bool s) { m_scale_to_mean = s; }

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
bool ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::GetScaleToMean() const { return m_scale_to_mean; }

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::SetPoolsize(const size_t n) { m_poolsize = n; }

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::SetSubregion(const TRegion &sr) { m_subregion = sr; m_hasSubregion = true; }

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::SetVerbose(const bool v) { m_verbose = v; }

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::SetOutputAllResiduals(const bool r) { m_allResiduals = r; }

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
RealTimeClock::TimeStampType ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::GetTotalTime() const { return m_elapsedTime; }

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
RealTimeClock::TimeStampType ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::GetMeanTime() const { return m_elapsedTime / m_unmaskedVoxels; }

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
SizeValueType ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::GetEvaluations() const { return m_unmaskedVoxels; }

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::SetInput(const size_t i, const TDataVectorImage *image) {
	if (i < m_algorithm->numInputs()) {
		this->SetNthInput(i, const_cast<TDataVectorImage*>(image));
	} else {
        itkExceptionMacro("Requested input " << i << " does not exist (" << m_algorithm->numInputs() << " inputs)");
	}
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::SetConst(const size_t i, const TScalarImage *image) {
	if (i < m_algorithm->numConsts()) {
		this->SetNthInput(m_algorithm->numInputs() + 1 + i, const_cast<TScalarImage*>(image));
	} else {
        itkExceptionMacro("Requested const input " << i << " does not exist (" << m_algorithm->numConsts() << " const inputs)");
	}
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::SetMask(const TScalarImage *image) {
	this->SetNthInput(m_algorithm->numInputs(), const_cast<TScalarImage*>(image));
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
auto ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::GetInput(const size_t i) const -> typename TDataVectorImage::ConstPointer {
	if (i < m_algorithm->numInputs()) {
		return static_cast<const TDataVectorImage *> (this->ProcessObject::GetInput(i));
	} else {
        itkExceptionMacro("Requested input " << i << " does not exist (" << m_algorithm->numInputs() << " inputs)");
	}
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
auto ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::GetConst(const size_t i) const -> typename TScalarImage::ConstPointer {
	if (i < m_algorithm->numConsts()) {
		size_t index = m_algorithm->numInputs() + 1 + i;
		return static_cast<const TScalarImage *> (this->ProcessObject::GetInput(index));
	} else {
        itkExceptionMacro("Requested const input " << i << " does not exist (" << m_algorithm->numConsts() << " const inputs)");
	}
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
auto ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::GetMask() const -> typename TScalarImage::ConstPointer {
	return static_cast<const TScalarImage *>(this->ProcessObject::GetInput(m_algorithm->numInputs()));
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
DataObject::Pointer ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::MakeOutput(unsigned int idx) {
    DataObject::Pointer output;
    if (idx == AllResidualsOutput) {
        auto img = TScalarVectorImage::New();
        output = img;
    } else if (idx == ResidualOutput) {
        auto img = TScalarImage::New();
        output = img;
    } else if (idx == IterationsOutput) {
        auto img = TIterationsImage::New();
        output = img;
    } else if (idx < (m_algorithm->numOutputs() + StartOutputs)) {
        output = (TScalarImage::New()).GetPointer();
    } else {
        itkExceptionMacro("Attempted to create output " << idx << ", this algorithm only has " << m_algorithm->numOutputs() << "+" << StartOutputs << " outputs.");
    }
    return output.GetPointer();
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
auto ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::GetOutput(const size_t i) -> TScalarImage *{
	if (i < m_algorithm->numOutputs()) {
        return dynamic_cast<TScalarImage *>(this->ProcessObject::GetOutput(i+StartOutputs));
	} else {
        itkExceptionMacro("Requested output " << std::to_string(i) << " is past maximum (" << std::to_string(m_algorithm->numOutputs()) << ")");
	}
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
auto ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::GetAllResidualsOutput() -> TScalarVectorImage *{
    return dynamic_cast<TScalarVectorImage *>(this->ProcessObject::GetOutput(AllResidualsOutput));
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
auto ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::GetResidualOutput() -> TScalarImage *{
    return dynamic_cast<TScalarImage *>(this->ProcessObject::GetOutput(ResidualOutput));
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
auto ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::GetIterationsOutput() -> TIterationsImage *{
    return dynamic_cast<TIterationsImage *>(this->ProcessObject::GetOutput(IterationsOutput));
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::GenerateOutputInformation() {
	Superclass::GenerateOutputInformation();
	size_t size = 0;
	for (size_t i = 0; i < m_algorithm->numInputs(); i++) {
		size += this->GetInput(i)->GetNumberOfComponentsPerPixel();
    }
	if (m_algorithm->dataSize() != size) {
        itkExceptionMacro("Sequence size (" << m_algorithm->dataSize() << ") does not match input size (" << size << ")");
	}
    if (size == 0) {
        itkExceptionMacro("Total input size cannot be 0");
    }

    auto input     = this->GetInput(0);
    auto region    = input->GetLargestPossibleRegion();
    auto spacing   = input->GetSpacing();
    auto origin    = input->GetOrigin();
    auto direction = input->GetDirection();
    if (m_verbose) std::cout << "Allocating output memory" << std::endl;
    for (size_t i = 0; i < m_algorithm->numOutputs(); i++) {
        auto op = this->GetOutput(i);
        op->SetRegions(region);
        op->SetSpacing(spacing);
        op->SetOrigin(origin);
        op->SetDirection(direction);
        op->Allocate(true);
    }
    if (m_allResiduals) {
        if (m_verbose) std::cout << "Allocating residuals memory" << std::endl;
        auto r = this->GetAllResidualsOutput();
        r->SetRegions(region);
        r->SetSpacing(spacing);
        r->SetOrigin(origin);
        r->SetDirection(direction);
        r->SetNumberOfComponentsPerPixel(size);
        r->Allocate(true);
    }
    auto r = this->GetResidualOutput();
    r->SetRegions(region);
    r->SetSpacing(spacing);
    r->SetOrigin(origin);
    r->SetDirection(direction);
    r->Allocate(true);
    auto i = this->GetIterationsOutput();
    i->SetRegions(region);
    i->SetSpacing(spacing);
    i->SetOrigin(origin);
    i->SetDirection(direction);
    i->Allocate(true);
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>::GenerateData() {
    const unsigned int LastDim = ImageDim - 1;
    TRegion region = this->GetInput(0)->GetLargestPossibleRegion();
    if (m_hasSubregion) {
        if (region.IsInside(m_subregion)) {
            region = m_subregion;
        } else {
            itkExceptionMacro("Specified subregion is not entirely inside image.");
        }
    }
    ImageRegionConstIterator<TScalarImage> maskIter;
    const auto mask = this->GetMask();
    if (mask) {
        if (m_verbose) std::cout << "Counting voxels in mask..." << std::endl;
        m_unmaskedVoxels = 0;
        maskIter = ImageRegionConstIterator<TScalarImage>(mask, region);
        maskIter.GoToBegin();
        while (!maskIter.IsAtEnd()) {
            if (maskIter.Get())
                ++m_unmaskedVoxels;
            ++maskIter;
        }
        maskIter.GoToBegin(); // Reset
        if (m_verbose) std::cout << "Found " << m_unmaskedVoxels << " unmasked voxels." << std::endl;
    } else {
        m_unmaskedVoxels = region.GetNumberOfPixels();
    }
	ProgressReporter progress(this, 0, m_unmaskedVoxels, 10);

    std::vector<ImageRegionConstIterator<TDataVectorImage>> dataIters(m_algorithm->numInputs());
	for (size_t i = 0; i < m_algorithm->numInputs(); i++) {
		dataIters[i] = ImageRegionConstIterator<TDataVectorImage>(this->GetInput(i), region);
	}

    std::vector<ImageRegionConstIterator<TScalarImage>> constIters(m_algorithm->numConsts());
	for (size_t i = 0; i < m_algorithm->numConsts(); i++) {
		typename TScalarImage::ConstPointer c = this->GetConst(i);
		if (c) {
			constIters[i] = ImageRegionConstIterator<TScalarImage>(c, region);
		}
	}
    std::vector<ImageRegionIterator<TScalarImage>> outputIters(m_algorithm->numOutputs());
	for (size_t i = 0; i < m_algorithm->numOutputs(); i++) {
		outputIters[i] = ImageRegionIterator<TScalarImage>(this->GetOutput(i), region);
	}
    ImageRegionIterator<TScalarVectorImage> allResidualsIter;
    if (m_allResiduals) {
        allResidualsIter = ImageRegionIterator<TScalarVectorImage>(this->GetAllResidualsOutput(), region);
    }
    ImageRegionIterator<TScalarImage> residualIter(this->GetResidualOutput(), region);
    ImageRegionIterator<TIterationsImage> iterationsIter(this->GetIterationsOutput(), region);
	typedef typename TAlgorithm::TArray TArray;
    if (m_verbose) std::cout << "Starting processing" << std::endl;
    QI::ThreadPool threadPool(m_poolsize);
    TimeProbe clock;
    clock.Start();
	while(!dataIters[0].IsAtEnd()) {
		if (!mask || maskIter.Get()) {
            auto task = [=] {
                TArray outputs = TArray::Zero(m_algorithm->numOutputs());
                TArray resids =  TArray::Zero(m_algorithm->dataSize());
                typename TAlgorithm::TIterations iterations{0};
                
                TArray constants = m_algorithm->defaultConsts();
                for (size_t i = 0; i < constIters.size(); i++) {
                    if (this->GetConst(i)) {
                        constants[i] = constIters[i].Get();
                    }
                }
                typename TAlgorithm::TInput allData(m_algorithm->dataSize());
                size_t dataIndex = 0;
                for (size_t i = 0; i < m_algorithm->numInputs(); i++) {
                    VariableLengthVector<TData> dataVector = dataIters[i].Get();
                    Eigen::Map<const Eigen::Array<TData, Eigen::Dynamic, 1>> data(dataVector.GetDataPointer(), dataVector.Size());
                    Eigen::Array<TData, Eigen::Dynamic, 1> scaled = data;
                    if (m_scale_to_mean) {
                        scaled /= scaled.mean();
                    }
                    allData.segment(dataIndex, data.rows()) = scaled.template cast<typename TAlgorithm::TScalar>();
                    dataIndex += data.rows();
                }
                m_algorithm->apply(allData, constants, outputs, resids, iterations);
                for (size_t i = 0; i < m_algorithm->numOutputs(); i++) {
                    outputIters[i].Set(static_cast<float>(outputs[i]));
                }
                Eigen::ArrayXf residF = resids.template cast<float>();
                VariableLengthVector<float> residVector(residF.data(), m_algorithm->dataSize());
                if (m_allResiduals) {
                    allResidualsIter.Set(residVector);
                }
                residualIter.Set(sqrt(resids.square().sum() / m_algorithm->dataSize()));
                iterationsIter.Set(iterations);
            };
            threadPool.enqueue(task);
            progress.CompletedPixel(); // We can get away with this because enqueue blocks if the queue is full
		} else {
            for (size_t i = 0; i < m_algorithm->numOutputs(); i++) {
                outputIters[i].Set(0);
            }
            VariableLengthVector<float> residZeros(m_algorithm->dataSize()); residZeros.Fill(0.);
            allResidualsIter.Set(residZeros);
            residualIter.Set(0);
            iterationsIter.Set(0);
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
			++outputIters[i];
		}
        if (m_allResiduals)
            ++allResidualsIter;
        ++residualIter;
        ++iterationsIter;
    }
    clock.Stop();
    m_elapsedTime = clock.GetTotal();
}
} // namespace ITK

#endif // APPLYALGORITHMFILTER_HXX
