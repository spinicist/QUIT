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

/*
template<typename TVImage, typename TAlgo>
void ApplyAlgorithmSliceBySliceFilter<TVImage, TAlgo>::GenerateInputRequestedRegion() {
	// call the superclass's implementation of this method, which
	// propagates the output requested region to all inputs
	Superclass::GenerateInputRequestedRegion();

	const TRegion &requestedInputRegion = this->GetDataInput(0)->GetRequestedRegion();
	// The requested region is the largest in all but the slice
	// dimension. In that dimension we can stream the requested
	// slices.
	TRegion inputRegion = this->GetDataInput(0)->GetLargestPossibleRegion();
	inputRegion.SetIndex(SliceDimension, requestedInputRegion.GetIndex(SliceDimension));
	inputRegion.SetSize(SliceDimension, requestedInputRegion.GetSize(SliceDimension));

	// Use the same requested region for each input, if an input image
	// is a different size and can't fulfill the request,
	// DataObject::PropagateRequestedRegion with throw
	for (int i = 0; i < this->m_algorithm->numInputs(); i++) {
		(this->GetDataInput(i))->SetRequestedRegion(inputRegion);
	}
	for (int i = 0; i < this->m_algorithm->numConsts(); i++) {
		if (this->GetConstInput(i))
			this->GetConstInput(i)->SetRequestedRegion(inputRegion);
	}
	if (this->GetMask())
		this->GetMask()->SetRequestedRegion(inputRegion);
}*/

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

	ProgressReporter progress(this, 0, requestedSize[SliceDimension]);

	// Allocate storage for slices
	//std::cout << "Setting storage for inputs" << std::endl;
	std::vector<typename TVectorSlice::Pointer> inputSlices(numInputs);
	for (unsigned int i = 0; i < numInputs; i++) {
		auto ip = this->GetDataInput(i);
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
		if (this->GetConstInput(i)) {
			const typename TImage::IndexType originIndex = this->GetConstInput(i)->GetRequestedRegion().GetIndex();
			typename TImage::PointType inputOrigin;
			this->GetDataInput(i)->TransformIndexToPhysicalPoint(originIndex, inputOrigin);

			typename TSlice::SpacingType sliceSpacing;
			typename TSlice::PointType sliceOrigin;
			unsigned int s = 0;
			for (unsigned int dim = 0; s < SliceDimension; ++dim, ++s) {
				if (dim == SliceDimension) {
					++dim;
				}
				sliceSpacing[s] = this->GetConstInput(i)->GetSpacing()[dim];
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
			const typename TImage::IndexType originIndex = this->GetMask()->GetRequestedRegion().GetIndex();
			typename TImage::PointType inputOrigin;
			this->GetMask()->TransformIndexToPhysicalPoint(originIndex, inputOrigin);

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

	//std::cout << "Starting" << std::endl;
	const int sliceRangeMax = (requestedSize[SliceDimension] + requestedIndex[SliceDimension]);
	for (m_sliceIndex = requestedIndex[SliceDimension]; m_sliceIndex < sliceRangeMax; ++m_sliceIndex ) {
		// say to the user that we are begining a new slice
		//std::cout << "Processing slice " << m_sliceIndex << std::endl;
		this->InvokeEvent(IterationEvent());

		typename TVectorImage::RegionType vectorInputRegion = this->GetInput(0)->GetRequestedRegion();
		typename TImage::RegionType inputRegion;
		vectorInputRegion.SetIndex(SliceDimension, m_sliceIndex);
		vectorInputRegion.SetSize(SliceDimension, 1);
		inputRegion.SetIndex(vectorInputRegion.GetIndex());
		inputRegion.SetSize(vectorInputRegion.GetSize());

		// this region is the current output region we are iterating on
		typename TVectorImage::RegionType vectorOutputRegion = this->GetOutput(0)->GetRequestedRegion();
		typename TImage::RegionType outputRegion;
		vectorOutputRegion.SetIndex(SliceDimension, m_sliceIndex);
		vectorOutputRegion.SetSize(SliceDimension, 1);
		outputRegion.SetIndex(vectorOutputRegion.GetIndex());
		outputRegion.SetSize(vectorOutputRegion.GetSize());

		//std::cout << "vectorInputRegion: " << vectorInputRegion << std::endl;
		//std::cout << "vectorSliceInputRegion: " << vectorSliceInputRegion << std::endl;
		//std::cout << "inputRegion: " << inputRegion << std::endl;
		//std::cout << "inputSliceRegion: " << sliceInputRegion << std::endl;

		//std::cout << "vectorOutputRegion: " << vectorOutputRegion << std::endl;
		//std::cout << "vectorSliceOutputRegion: " << vectorSliceOutputRegion << std::endl;
		//std::cout << "outputRegion: " << outputRegion << std::endl;
		//std::cout << "outputSliceRegion: " << sliceOutputRegion << std::endl;

		itkAssertOrThrowMacro( vectorInputRegion.GetNumberOfPixels() == vectorSliceInputRegion.GetNumberOfPixels(), "inputRegion.GetNumberOfPixels() == internalInputRegion.GetNumberOfPixel()" );
		itkAssertOrThrowMacro( vectorOutputRegion.GetNumberOfPixels() == vectorSliceOutputRegion.GetNumberOfPixels(), "outputRegion.GetNumberOfPixels() == internalOutputRegion.GetNumberOfPixel()" );

		//std::cout << "Creating filter" << std::endl;
		typename TSliceFilter::Pointer sliceFilter = TSliceFilter::New();
		sliceFilter->SetAlgorithm(this->m_algorithm);

		// reallocate the internal input at each slice, so the slice by slice filter can work
		// even if the pipeline is run in place
		//std::cout << "Setting input regions" << std::endl;
		for (int i = 0; i < numInputs; i++) {
			inputSlices[i]->SetRegions(vectorSliceInputRegion);
			inputSlices[i]->Allocate();
			sliceFilter->SetDataInput(i, inputSlices[i]);
			ImageAlgorithm::Copy(this->GetDataInput(i).GetPointer(), inputSlices[i].GetPointer(), vectorInputRegion, vectorSliceInputRegion);
		}
		//std::cout << "Setting Const regions" << std::endl;
		for (int i = 0; i < numConsts; i++) {
			if (this->GetConstInput(i)) {
				constSlices[i]->SetRegions(sliceInputRegion);
				constSlices[i]->Allocate();
				sliceFilter->SetConstInput(i, constSlices[i]);
				ImageAlgorithm::Copy(this->GetConstInput(i).GetPointer(), constSlices[i].GetPointer(), inputRegion, sliceInputRegion);
			}
		}
		//std::cout << "Setting mask region" << std::endl;
		if (this->GetMask()) {
			maskSlice->SetRegions(sliceInputRegion);
			maskSlice->Allocate();
			sliceFilter->SetMask(maskSlice);
			ImageAlgorithm::Copy(this->GetMask().GetPointer(), maskSlice.GetPointer(), inputRegion, sliceInputRegion);
		}

		/*for (int i = 0; i < numOutputs; i++) {
			sliceFilter->GetOutput(i)->SetRequestedRegion(vectorSliceOutputRegion);
		}
		sliceFilter->GetResidOutput()->SetRequestedRegion(vectorSliceOutputRegion);*/

		//std::cout << "Updating" << std::endl;
		sliceFilter->Update();
		//std::cout << "Finished updating" << std::endl;
		// and copy the output slice to the output image
		for (int i = 0; i < numOutputs; i++) {
			//std::cout << "Copying output " << i << std::endl;
			ImageAlgorithm::Copy(sliceFilter->GetOutput(i), this->GetOutput(i), vectorSliceOutputRegion, vectorOutputRegion);
		}
		//std::cout << "Copying resid" << std::endl;
		ImageAlgorithm::Copy(sliceFilter->GetResidOutput(), this->GetResidOutput(), vectorSliceOutputRegion, vectorOutputRegion);
		//std::cout << "Completing pixel" << std::endl;
		progress.CompletedPixel();
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
