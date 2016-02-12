#ifndef IMAGETOVECTORFILTER_HXX
#define IMAGETOVECTORFILTER_HXX

namespace itk {

template<typename TInput>
ImageToVectorFilter<TInput>::ImageToVectorFilter() {
	m_compose = ComposeType::New();
    m_BlockStart = 0;
    m_BlockSize = 0;
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
    if (m_BlockSize == 0) {
        m_BlockSize = inputRegion.GetSize()[OutputDimension];
    } else if (m_BlockSize > inputRegion.GetSize()[OutputDimension]) {
        itkExceptionMacro("Block size is larger than input image length.");
    } else if (inputRegion.GetSize()[OutputDimension] % m_BlockSize) {
        itkExceptionMacro("Block size does not divide input image length.");
    }
    outputPtr->SetNumberOfComponentsPerPixel(m_BlockSize);
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
    //std::cout << "END " << __PRETTY_FUNCTION__ << std::endl;
}

template<typename TInput>
void ImageToVectorFilter<TInput>::GenerateData() {
    //std::cout << __PRETTY_FUNCTION__ << std::endl;
	auto input = this->GetInput();
	auto region = input->GetLargestPossibleRegion();
    size_t blockEnd = m_BlockStart + m_BlockSize;
    size_t inputLength = region.GetSize()[OutputDimension];
    if (blockEnd > inputLength) {
        itkExceptionMacro("Block end " << blockEnd << " would be greater than input length (" << inputLength << ")");
    }
    region.GetModifiableSize()[OutputDimension] = 0;
    for (int i = 0; i < m_BlockSize; i ++) {
        region.GetModifiableIndex()[OutputDimension] = i + m_BlockStart;
        //std::cout << region << std::endl;
		auto volume = ExtractType::New();
		volume->SetExtractionRegion(region);
		volume->SetInput(input);
		volume->SetDirectionCollapseToSubmatrix();
		volume->Update();
		m_compose->SetInput(i, volume->GetOutput());
	}
	m_compose->Update();
	this->GraftOutput(m_compose->GetOutput());
    //std::cout << "END " << __PRETTY_FUNCTION__ << std::endl;
}

} // End namespace itk

#endif // Define IMAGETOVECTORFILTER_HXX
