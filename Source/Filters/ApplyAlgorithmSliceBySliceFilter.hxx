#ifndef APPLYALGORITHMSLICEBYSLICEFILTER_HXX
#define APPLYALGORITHMSLICEBYSLICEFILTER_HXX

namespace itk {

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
ApplyAlgorithmSliceBySliceFilter<TAlgorithm, TData, TScalar, ImageDim>::ApplyAlgorithmSliceBySliceFilter() :
	ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim>() {
	//std::cout << __PRETTY_FUNCTION__ << std::endl;
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmSliceBySliceFilter<TAlgorithm, TData, TScalar, ImageDim>::VerifyInputInformation() {
	//std::cout << __PRETTY_FUNCTION__ << std::endl;
	Superclass::VerifyInputInformation();

	if (!this->m_algorithm) {
		itkExceptionMacro("Algorithm must be set.");
	}
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmSliceBySliceFilter<TAlgorithm, TData, TScalar, ImageDim>::SetSlices(const int start, const int stop) {
	m_startSlice = start;
	m_stopSlice = stop;
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
int ApplyAlgorithmSliceBySliceFilter<TAlgorithm, TData, TScalar, ImageDim>::GetSliceIndex() const { return m_sliceIndex; }

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmSliceBySliceFilter<TAlgorithm, TData, TScalar, ImageDim>::GenerateData() {
	//std::cout << __PRETTY_FUNCTION__ << std::endl;
	const int numInputs = this->m_algorithm->numInputs();
	const int numConsts = this->m_algorithm->numConsts();
	const int numOutputs = this->m_algorithm->numOutputs();

	//this->AllocateOutputs();

	const typename TDataVectorImage::RegionType requestedRegion = this->GetOutput( 0 )->GetRequestedRegion();
	const typename TDataVectorImage::RegionType::IndexType requestedIndex = requestedRegion.GetIndex();
	const typename TDataVectorImage::RegionType::SizeType requestedSize = requestedRegion.GetSize();

	typename TDataVectorSlice::RegionType vectorSliceInputRegion, vectorSliceOutputRegion;
	typename TScalarSlice::RegionType sliceInputRegion, sliceOutputRegion;
	// copy the requested region to the internal slice region in dimension order
	unsigned int slice_i = 0;
	for (unsigned int i = 0; slice_i < SliceDim; ++i, ++slice_i) {
		if (i == SliceDim ) {
			++i;
		}
		vectorSliceOutputRegion.SetSize(slice_i, requestedSize[i]);
		vectorSliceOutputRegion.SetIndex(slice_i, requestedIndex[i]);
		vectorSliceInputRegion.SetSize(slice_i, this->GetInput(0)->GetRequestedRegion().GetSize(i));
		vectorSliceInputRegion.SetIndex(slice_i, this->GetInput(0)->GetRequestedRegion().GetIndex(i));
	}
	sliceInputRegion.SetSize(vectorSliceInputRegion.GetSize());
	sliceInputRegion.SetIndex(vectorSliceInputRegion.GetIndex());
	sliceOutputRegion.SetSize(vectorSliceOutputRegion.GetSize());
	sliceOutputRegion.SetIndex(vectorSliceOutputRegion.GetIndex());

	// Allocate storage for slices
	//std::cout << "Setting storage for inputs" << std::endl;
	std::vector<typename TDataVectorSlice::Pointer> inputSlices(numInputs);
	for (unsigned int i = 0; i < numInputs; i++) {
		const auto ip = this->GetInput(i);
		const typename TDataVectorImage::IndexType originIndex = ip->GetRequestedRegion().GetIndex();
		typename TDataVectorImage::PointType inputOrigin;
		ip->TransformIndexToPhysicalPoint(originIndex, inputOrigin);

		typename TDataVectorSlice::SpacingType sliceSpacing;
		typename TDataVectorSlice::PointType sliceOrigin;
		unsigned int s = 0;
		for (unsigned int dim = 0; s < SliceDim; ++dim, ++s) {
			if (dim == SliceDim) {
				++dim;
			}
			sliceSpacing[s] = ip->GetSpacing()[dim];
			sliceOrigin[s] = inputOrigin[dim];
		}
		inputSlices[i] = TDataVectorSlice::New();
		inputSlices[i]->SetSpacing(sliceSpacing);
		inputSlices[i]->SetOrigin(sliceOrigin);
		inputSlices[i]->SetNumberOfComponentsPerPixel(ip->GetNumberOfComponentsPerPixel());
	}

	//std::cout << "Setting storage for consts" << std::endl;
	std::vector<typename TScalarSlice::Pointer> constSlices(numConsts);
	for (unsigned int i = 0; i < numConsts; i++) {
		const auto cp = this->GetConst(i);
		if (cp) {
			const typename TScalarImage::IndexType originIndex = cp->GetRequestedRegion().GetIndex();
			typename TScalarImage::PointType inputOrigin;
			cp->TransformIndexToPhysicalPoint(originIndex, inputOrigin);

			typename TScalarSlice::SpacingType sliceSpacing;
			typename TScalarSlice::PointType sliceOrigin;
			unsigned int s = 0;
			for (unsigned int dim = 0; s < SliceDim; ++dim, ++s) {
				if (dim == SliceDim) {
					++dim;
				}
				sliceSpacing[s] = this->GetConst(i)->GetSpacing()[dim];
				sliceOrigin[s] = inputOrigin[dim];
			}
			constSlices[i] = TScalarSlice::New();
			constSlices[i]->SetSpacing(sliceSpacing);
			constSlices[i]->SetOrigin(sliceOrigin);
		} else {
			constSlices[i] = ITK_NULLPTR;
		}
	}

	//std::cout << "Setting storage for mask" << std::endl;
	typename TScalarSlice::Pointer maskSlice = ITK_NULLPTR;
	if (this->GetMask()) {
		const auto mp = this->GetMask();
		const typename TScalarImage::IndexType originIndex = mp->GetRequestedRegion().GetIndex();
		typename TScalarImage::PointType inputOrigin;
		mp->TransformIndexToPhysicalPoint(originIndex, inputOrigin);

		typename TScalarSlice::SpacingType sliceSpacing;
		typename TScalarSlice::PointType sliceOrigin;
		unsigned int s = 0;
		for (unsigned int dim = 0; s < SliceDim; ++dim, ++s) {
			if (dim == SliceDim) {
				++dim;
			}
			sliceSpacing[s] = this->GetMask()->GetSpacing()[dim];
			sliceOrigin[s] = inputOrigin[dim];
		}
		maskSlice = TScalarSlice::New();
		maskSlice->SetSpacing(sliceSpacing);
		maskSlice->SetOrigin(sliceOrigin);
	}

	for (int i = 0; i < numOutputs; i++) {
		this->GetOutput(i)->FillBuffer(0);
	}
	VariableLengthVector<TScalar> zero(this->GetResidOutput()->GetNumberOfComponentsPerPixel());
	zero.Fill(0);
	this->GetResidOutput()->FillBuffer(zero);

	//std::cout << "Starting" << std::endl;
	const int sliceRangeMax = (requestedSize[SliceDim] + requestedIndex[SliceDim]);
	if (m_stopSlice == 0)
		m_stopSlice = sliceRangeMax;
	ProgressReporter progress(this, 0, sliceRangeMax);
	for (m_sliceIndex = requestedIndex[SliceDim]; m_sliceIndex < sliceRangeMax; ++m_sliceIndex ) {
		// this region is the current output region we are iterating on
		typename TScalarVectorImage::RegionType vectorOutputRegion = this->GetOutput(0)->GetRequestedRegion();
		typename TScalarImage::RegionType outputRegion;
		vectorOutputRegion.SetIndex(SliceDim, m_sliceIndex);
		vectorOutputRegion.SetSize(SliceDim, 1);
		outputRegion.SetIndex(vectorOutputRegion.GetIndex());
		outputRegion.SetSize(vectorOutputRegion.GetSize());

		if ((m_sliceIndex >= m_startSlice) && (m_sliceIndex < m_stopSlice)) {
			this->InvokeEvent(IterationEvent());
			typename TDataVectorImage::RegionType vectorInputRegion = this->GetInput(0)->GetRequestedRegion();
			typename TScalarImage::RegionType inputRegion;
			vectorInputRegion.SetIndex(SliceDim, m_sliceIndex);
			vectorInputRegion.SetSize(SliceDim, 1);
			inputRegion.SetIndex(vectorInputRegion.GetIndex());
			inputRegion.SetSize(vectorInputRegion.GetSize());

			itkAssertOrThrowMacro(vectorInputRegion.GetNumberOfPixels() == vectorSliceInputRegion.GetNumberOfPixels(), "inputRegion.GetNumberOfPixels() == internalInputRegion.GetNumberOfPixel()");
			itkAssertOrThrowMacro(vectorOutputRegion.GetNumberOfPixels() == vectorSliceOutputRegion.GetNumberOfPixels(), "outputRegion.GetNumberOfPixels() == internalOutputRegion.GetNumberOfPixel()");

			//std::cout << "Creating filter" << std::endl;
			typename TSliceFilter::Pointer sliceFilter = TSliceFilter::New();
			sliceFilter->SetAlgorithm(this->m_algorithm);
			sliceFilter->SetScaleToMean(this->m_scale_to_mean);

			// reallocate the internal input at each slice, so the slice by slice filter can work
			// even if the pipeline is run in place
			//std::cout << "Setting input regions" << std::endl;
			for (int i = 0; i < numInputs; i++) {
				inputSlices[i]->SetRegions(vectorSliceInputRegion);
				inputSlices[i]->Allocate();
				sliceFilter->SetInput(i, inputSlices[i]);
				ImageAlgorithm::Copy(this->GetInput(i).GetPointer(), inputSlices[i].GetPointer(), vectorInputRegion, vectorSliceInputRegion);
			}
			//std::cout << "Setting Const regions" << std::endl;
			for (int i = 0; i < numConsts; i++) {
				if (this->GetConst(i)) {
					constSlices[i]->SetRegions(sliceInputRegion);
					constSlices[i]->Allocate();
					sliceFilter->SetConst(i, constSlices[i]);
					ImageAlgorithm::Copy(this->GetConst(i).GetPointer(), constSlices[i].GetPointer(), inputRegion, sliceInputRegion);
				}
			}
			//std::cout << "Setting mask region" << std::endl;
			if (this->GetMask()) {
				maskSlice->SetRegions(sliceInputRegion);
				maskSlice->Allocate();
				sliceFilter->SetMask(maskSlice);
				ImageAlgorithm::Copy(this->GetMask().GetPointer(), maskSlice.GetPointer(), inputRegion, sliceInputRegion);
			}

			sliceFilter->Update();

			// and copy the output slice to the output image
			for (int i = 0; i < numOutputs; i++) {
				ImageAlgorithm::Copy(sliceFilter->GetOutput(i), this->GetOutput(i), sliceOutputRegion, outputRegion);
			}
			ImageAlgorithm::Copy(sliceFilter->GetResidOutput(), this->GetResidOutput(), vectorSliceOutputRegion, vectorOutputRegion);
			progress.CompletedPixel(); // Put this here so skipped slices don't get a message
		}
		//std::cout << "Finised loop" << std::endl;
	}
	//std::cout << "Finished" << std::endl;
}

template<typename TAlgorithm, typename TData, typename TScalar, unsigned int ImageDim>
void ApplyAlgorithmSliceBySliceFilter<TAlgorithm, TData, TScalar, ImageDim>::PrintSelf(std::ostream & os, Indent indent) const {
	Superclass::PrintSelf(os, indent);

	os << indent << "Dimension: " << SliceDim << std::endl;
	os << indent << "SliceIndex: " << m_sliceIndex << std::endl;
}

} // namespace ITK

#endif // APPLYALGORITHMSLICEBYSLICEFILTER_HXX
