template<typename TPixel>
VectorToImageFilter<TPixel>::VectorToImageFilter() {
	//std::cout << __PRETTY_FUNCTION__ << std::endl;
	m_tiler = TTiler::New();
}

template<typename TPixel>
void VectorToImageFilter<TPixel>::GenerateOutputInformation() {
	//std::cout << __PRETTY_FUNCTION__ << std::endl;
	auto out = this->GetOutput();
	auto in  = this->GetInput();
	if ( !out || !in ) {
		return;
	}

	itk::FixedArray<unsigned int, 4> layout;
	layout[0] = layout[1] = layout[2] = 1;
	layout[3] = in->GetNumberOfComponentsPerPixel();
	m_tiler->SetLayout(layout);
	m_indexers.resize(in->GetNumberOfComponentsPerPixel());
	for (size_t i = 0; i < in->GetNumberOfComponentsPerPixel(); i++) {
		m_indexers[i] = TIndexer::New();
		m_indexers[i]->SetIndex(i);
		m_indexers[i]->SetInput(in);
		m_tiler->SetInput(i, m_indexers[i]->GetOutput());
	}
	m_tiler->UpdateLargestPossibleRegion();
	this->GraftOutput(m_tiler->GetOutput());
}

template<typename TPixel>
void VectorToImageFilter<TPixel>::Update() {
	//std::cout << __PRETTY_FUNCTION__ << std::endl;
	m_tiler->Update();
	//this->GraftOutput(m_tiler->GetOutput());
}
