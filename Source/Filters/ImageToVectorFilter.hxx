#ifndef IMAGETOVECTORFILTER_HXX
#define IMAGETOVECTORFILTER_HXX

namespace itk {

template<typename TInput>
ImageToVectorFilter<TInput>::ImageToVectorFilter() {
	m_compose = ComposeType::New();
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
	outputPtr->SetNumberOfComponentsPerPixel(inputPtr->GetLargestPossibleRegion().GetSize()[OutputDimension]);
	//std::cout << "End " << __PRETTY_FUNCTION__ << std::endl;
}

template<typename TInput>
void ImageToVectorFilter<TInput>::GenerateData() {
	auto input = this->GetInput();
	auto region = input->GetLargestPossibleRegion();
	int size = region.GetSize()[OutputDimension];
	region.GetModifiableSize()[OutputDimension] = 0;
	for (int i = 0; i < size; i ++) {
		region.GetModifiableIndex()[OutputDimension] = i;
		auto volume = ExtractType::New();
		volume->SetExtractionRegion(region);
		volume->SetInput(input);
		volume->SetDirectionCollapseToSubmatrix();
		volume->Update();
		m_compose->SetInput(i, volume->GetOutput());
	}
	m_compose->Update();
	this->GraftOutput(m_compose->GetOutput());
}

} // End namespace itk

#endif // Define IMAGETOVECTORFILTER_HXX
