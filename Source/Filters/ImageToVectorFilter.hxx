#include "ImageToVectorFilter.h"

template<typename PixelType>
ImageToVectorFilter<PixelType>::ImageToVectorFilter() {
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	m_compose = ComposeType::New();
	//this->GraftOutput(m_compose->GetOutput());
}

template<typename PixelType>
void ImageToVectorFilter<PixelType>::GenerateOutputInformation() {
	std::cout << __PRETTY_FUNCTION__ << std::endl;
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
	std::cout << "inputPtr->GetDirection()" << std::endl << inputPtr->GetDirection() << std::endl;
	std::cout << "outputPtr->GetDirection()" << std::endl << outputPtr->GetDirection() << std::endl;
}

template<typename PixelType>
void ImageToVectorFilter<PixelType>::GenerateData() {
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	auto input = this->GetInput();
	auto region = input->GetLargestPossibleRegion();
	size_t nVols = region.GetSize()[3];
	region.GetModifiableSize()[3] = 0;
	std::cout << "input" << std::endl << input->GetDirection() << std::endl;
	for (int i = 0; i < nVols; i++) {
		region.GetModifiableIndex()[3] = i;
		auto volume = ExtractType::New();
		volume->SetExtractionRegion(region);
		volume->SetInput(input);
		std::cout << "volume before " << i << std::endl << volume->GetOutput()->GetDirection() << std::endl;
		volume->SetDirectionCollapseToSubmatrix();
		std::cout << "volume after" << i << std::endl << volume->GetOutput()->GetDirection() << std::endl;
		volume->Update();
		m_compose->SetInput(i, volume->GetOutput());
	}
	m_compose->Update();
	std::cout << "m_compose" << std::endl << m_compose->GetOutput()->GetDirection() << std::endl;
	std::cout << "output before" << std::endl << this->GetOutput()->GetDirection() << std::endl;
	this->GraftOutput(m_compose->GetOutput());
	std::cout << "output after" << std::endl << this->GetOutput()->GetDirection() << std::endl;
}
