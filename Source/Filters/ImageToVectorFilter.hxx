template<typename PixelType>
ImageToVectorFilter<PixelType>::ImageToVectorFilter() {
	m_compose = ComposeType::New();
}

template<typename PixelType>
void ImageToVectorFilter<PixelType>::GenerateOutputInformation() {
	typename Superclass::OutputImagePointer outputPtr = this->GetOutput();
	typename Superclass::InputImageConstPointer inputPtr  = this->GetInput();
	if ( !outputPtr || !inputPtr ) {
		return;
	}
	typename InputType::RegionType inputRegion = inputPtr->GetLargestPossibleRegion();
	outputPtr->SetLargestPossibleRegion(inputRegion.Slice(3));

	typename InputType::SpacingType spacing = inputPtr->GetSpacing();
	typename InputType::PointType   origin  = inputPtr->GetOrigin();
	typename InputType::DirectionType direction = inputPtr->GetDirection();
	typename OutputType::SpacingType outSpacing;
	typename OutputType::PointType   outOrigin;
	typename OutputType::DirectionType outDirection;

	for (int i = 0; i < 3; i++) {
		outSpacing[i] = spacing[i];
		outOrigin[i] =  origin[i];
		for (int j = 0; j < 3; j++) {
			outDirection[i][j] = direction[i][j];
		}
	}
	outputPtr->SetSpacing(outSpacing);
	outputPtr->SetOrigin(outOrigin);
	outputPtr->SetDirection(outDirection);
	outputPtr->SetNumberOfComponentsPerPixel(inputRegion.GetSize()[3]);
}

template<typename PixelType>
void ImageToVectorFilter<PixelType>::GenerateData() {
	auto input = this->GetInput();
	auto region = input->GetLargestPossibleRegion();
	size_t nVols = region.GetSize()[3];
	region.GetModifiableSize()[3] = 0;
	for (int i = 0; i < nVols; i++) {
		region.GetModifiableIndex()[3] = i;
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
