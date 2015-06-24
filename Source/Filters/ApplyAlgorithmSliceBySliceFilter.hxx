#ifndef APPLYALGORITHMSLICEBYSLICEFILTER_HXX
#define APPLYALGORITHMSLICEBYSLICEFILTER_HXX

namespace itk {

template<typename TVImage, typename TAlgo>
ApplyAlgorithmSliceBySliceFilter<TVImage, TAlgo>::ApplyAlgorithmSliceBySliceFilter() :
	ApplyAlgorithmFilter<TVImage, TAlgo>() {
	//std::cout << __PRETTY_FUNCTION__ << std::endl;
}

template<typename TVImage, typename TAlgo>
void ApplyAlgorithmSliceBySliceFilter<TVImage, TAlgo>::VerifyInputInformation() {
	//std::cout << __PRETTY_FUNCTION__ << std::endl;
	Superclass::VerifyInputInformation();

	if (!this->m_algorithm) {
		itkExceptionMacro("Algorithm must be set.");
	}
}

template<typename TVImage, typename TAlgo>
void ApplyAlgorithmSliceBySliceFilter<TVImage, TAlgo>::SetSlices(const int start, const int stop) {
	m_startSlice = start;
	m_stopSlice = stop;
}

template<typename TVImage, typename TAlgo>
int ApplyAlgorithmSliceBySliceFilter<TVImage, TAlgo>::GetSliceIndex() const { return m_sliceIndex; }

template<typename TVImage, typename TAlgo>
void ApplyAlgorithmSliceBySliceFilter<TVImage, TAlgo>::GenerateData() {
	//std::cout << __PRETTY_FUNCTION__ << std::endl;
	const int numInputs = this->m_algorithm->numInputs();
	const int numConsts = this->m_algorithm->numConsts();
	const int numOutputs = this->m_algorithm->numOutputs();

	//this->AllocateOutputs();

	const typename TVectorImage::RegionType requestedRegion = this->GetOutput( 0 )->GetRequestedRegion();
	const typename TVectorImage::RegionType::IndexType requestedIndex = requestedRegion.GetIndex();
	const typename TVectorImage::RegionType::SizeType requestedSize = requestedRegion.GetSize();

	typename TVectorSlice::RegionType vectorSliceInputRegion, vectorSliceOutputRegion;
	typename TSlice::RegionType sliceInputRegion, sliceOutputRegion;
	// copy the requested region to the internal slice region in dimension order
	unsigned int slice_i = 0;
	for (unsigned int i = 0; slice_i < SliceDimension; ++i, ++slice_i) {
		if (i == SliceDimension ) {
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
	std::vector<typename TVectorSlice::Pointer> inputSlices(numInputs);
	for (unsigned int i = 0; i < numInputs; i++) {
		const auto ip = this->GetInput(i);
		const typename TVectorImage::IndexType originIndex = ip->GetRequestedRegion().GetIndex();
		typename TVectorImage::PointType inputOrigin;
		ip->TransformIndexToPhysicalPoint(originIndex, inputOrigin);

		typename TVectorSlice::SpacingType sliceSpacing;
		typename TVectorSlice::PointType sliceOrigin;
		unsigned int s = 0;
		for (unsigned int dim = 0; s < SliceDimension; ++dim, ++s) {
			if (dim == SliceDimension) {
				++dim;
			}
			sliceSpacing[s] = ip->GetSpacing()[dim];
			sliceOrigin[s] = inputOrigin[dim];
		}
		inputSlices[i] = TVectorSlice::New();
		inputSlices[i]->SetSpacing(sliceSpacing);
		inputSlices[i]->SetOrigin(sliceOrigin);
		inputSlices[i]->SetNumberOfComponentsPerPixel(ip->GetNumberOfComponentsPerPixel());
	}

	//std::cout << "Setting storage for consts" << std::endl;
	std::vector<typename TSlice::Pointer> constSlices(numConsts);
	for (unsigned int i = 0; i < numConsts; i++) {
		const auto cp = this->GetConst(i);
		if (cp) {
			const typename TImage::IndexType originIndex = cp->GetRequestedRegion().GetIndex();
			typename TImage::PointType inputOrigin;
			cp->TransformIndexToPhysicalPoint(originIndex, inputOrigin);

			typename TSlice::SpacingType sliceSpacing;
			typename TSlice::PointType sliceOrigin;
			unsigned int s = 0;
			for (unsigned int dim = 0; s < SliceDimension; ++dim, ++s) {
				if (dim == SliceDimension) {
					++dim;
				}
				sliceSpacing[s] = this->GetConst(i)->GetSpacing()[dim];
				sliceOrigin[s] = inputOrigin[dim];
			}
			constSlices[i] = TSlice::New();
			constSlices[i]->SetSpacing(sliceSpacing);
			constSlices[i]->SetOrigin(sliceOrigin);
		} else {
			constSlices[i] = ITK_NULLPTR;
		}
	}

	//std::cout << "Setting storage for mask" << std::endl;
	typename TSlice::Pointer maskSlice = ITK_NULLPTR;
	if (this->GetMask()) {
		const auto mp = this->GetMask();
		const typename TImage::IndexType originIndex = mp->GetRequestedRegion().GetIndex();
		typename TImage::PointType inputOrigin;
		mp->TransformIndexToPhysicalPoint(originIndex, inputOrigin);

		typename TSlice::SpacingType sliceSpacing;
		typename TSlice::PointType sliceOrigin;
		unsigned int s = 0;
		for (unsigned int dim = 0; s < SliceDimension; ++dim, ++s) {
			if (dim == SliceDimension) {
				++dim;
			}
			sliceSpacing[s] = this->GetMask()->GetSpacing()[dim];
			sliceOrigin[s] = inputOrigin[dim];
		}
		maskSlice = TSlice::New();
		maskSlice->SetSpacing(sliceSpacing);
		maskSlice->SetOrigin(sliceOrigin);
	}

	for (int i = 0; i < numOutputs; i++) {
		this->GetOutput(i)->FillBuffer(0);
	}
	TVector zero(this->GetResidOutput()->GetNumberOfComponentsPerPixel());
	zero.Fill(0);
	this->GetResidOutput()->FillBuffer(zero);

	//std::cout << "Starting" << std::endl;
	const int sliceRangeMax = (requestedSize[SliceDimension] + requestedIndex[SliceDimension]);
	if (m_stopSlice == 0)
		m_stopSlice = sliceRangeMax;
	ProgressReporter progress(this, 0, sliceRangeMax);
	for (m_sliceIndex = requestedIndex[SliceDimension]; m_sliceIndex < sliceRangeMax; ++m_sliceIndex ) {
		// this region is the current output region we are iterating on
		typename TVectorImage::RegionType vectorOutputRegion = this->GetOutput(0)->GetRequestedRegion();
		typename TImage::RegionType outputRegion;
		vectorOutputRegion.SetIndex(SliceDimension, m_sliceIndex);
		vectorOutputRegion.SetSize(SliceDimension, 1);
		outputRegion.SetIndex(vectorOutputRegion.GetIndex());
		outputRegion.SetSize(vectorOutputRegion.GetSize());

		if ((m_sliceIndex >= m_startSlice) && (m_sliceIndex < m_stopSlice)) {
			this->InvokeEvent(IterationEvent());
			typename TVectorImage::RegionType vectorInputRegion = this->GetInput(0)->GetRequestedRegion();
			typename TImage::RegionType inputRegion;
			vectorInputRegion.SetIndex(SliceDimension, m_sliceIndex);
			vectorInputRegion.SetSize(SliceDimension, 1);
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
}

template< typename TImage, typename TAlgo>
void ApplyAlgorithmSliceBySliceFilter<TImage, TAlgo>::PrintSelf(std::ostream & os, Indent indent) const {
	Superclass::PrintSelf(os, indent);

	os << indent << "Dimension: " << SliceDimension << std::endl;
	os << indent << "SliceIndex: " << m_sliceIndex << std::endl;
}

} // namespace ITK

#endif // APPLYALGORITHMSLICEBYSLICEFILTER_HXX
