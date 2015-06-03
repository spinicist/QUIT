#ifndef IMAGETOVECTORFILTER_HXX
#define IMAGETOVECTORFILTER_HXX

namespace itk {

template<typename TInput>
ImageToVectorFilter<TInput>::ImageToVectorFilter() {
	m_compose = ComposeType::New();
}

template<typename TInput>
void ImageToVectorFilter<TInput>::SetStartStop(size_t start, size_t stop) {
	m_start = start;
	m_stop = stop;
}

template<typename TInput>
void ImageToVectorFilter<TInput>::SetStride(size_t s) {
	m_stride = s;
}

template<typename TInput>
void ImageToVectorFilter<TInput>::GenerateOutputInformation() {
	//std::cout << __PRETTY_FUNCTION__ << std::endl;
	typename Superclass::OutputImagePointer outputPtr = this->GetOutput();
	typename Superclass::InputImageConstPointer inputPtr  = this->GetInput();
	if ( !outputPtr || !inputPtr ) {
		return;
	}
	typename TInput::RegionType inputRegion = inputPtr->GetLargestPossibleRegion();
	outputPtr->SetLargestPossibleRegion(inputRegion.Slice(OutputDimension));

	typename TInput::SpacingType spacing = inputPtr->GetSpacing();
	typename TInput::PointType   origin  = inputPtr->GetOrigin();
	typename TInput::DirectionType direction = inputPtr->GetDirection();
	typename TOutput::SpacingType outSpacing;
	typename TOutput::PointType   outOrigin;
	typename TOutput::DirectionType outDirection;

	for (int i = 0; i < (OutputDimension); i++) {
		outSpacing[i] = spacing[i];
		outOrigin[i] =  origin[i];
		for (int j = 0; j < (OutputDimension); j++) {
			outDirection[i][j] = direction[i][j];
		}
	}
	outputPtr->SetSpacing(outSpacing);
	outputPtr->SetOrigin(outOrigin);
	outputPtr->SetDirection(outDirection);

	if (m_stop == 0)
		m_stop = inputPtr->GetLargestPossibleRegion().GetSize()[OutputDimension];
	m_size = (m_stop - m_start) / m_stride;
	outputPtr->SetNumberOfComponentsPerPixel(m_size);
	//std::cout << "End " << __PRETTY_FUNCTION__ << std::endl;
}

template<typename TInput>
void ImageToVectorFilter<TInput>::GenerateData() {
	auto input = this->GetInput();
	auto region = input->GetLargestPossibleRegion();
	region.GetModifiableSize()[OutputDimension] = 0;
	int outI = 0;
	for (int i = m_start; i < m_stop; i += m_stride) {
		region.GetModifiableIndex()[OutputDimension] = i;
		auto volume = ExtractType::New();
		volume->SetExtractionRegion(region);
		volume->SetInput(input);
		volume->SetDirectionCollapseToSubmatrix();
		volume->Update();
		m_compose->SetInput(outI++, volume->GetOutput());
	}
	m_compose->Update();
	this->GraftOutput(m_compose->GetOutput());
}

} // End namespace itk

#endif // Define IMAGETOVECTORFILTER_HXX
