#ifndef VECTORTOIMAGEFILTER_HXX
#define VECTORTOIMAGEFILTER_HXX

namespace itk {

template<typename TPixel>
VectorToImageFilter<TPixel>::VectorToImageFilter() {
	m_tiler = TTiler::New();
}

template<typename TPixel>
void VectorToImageFilter<TPixel>::GenerateOutputInformation() {
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
	auto spacing    = in->GetSpacing();
	auto origin     = in->GetOrigin();
	auto direction  = in->GetDirection();
	auto outSpacing   = out->GetSpacing();
	auto outOrigin    = out->GetOrigin();
	auto outDirection = out->GetDirection();
	outSpacing.Fill(1);
	outOrigin.Fill(1);
	outDirection.SetIdentity();
	for (int i = 0; i < (InputDimension); i++) {
		outSpacing[i] = spacing[i];
		outOrigin[i] =  origin[i];
		for (int j = 0; j < (InputDimension); j++) {
			outDirection[i][j] = direction[i][j];
		}
	}
	out->SetSpacing(outSpacing);
	out->SetOrigin(outOrigin);
	out->SetDirection(outDirection);
}

template<typename TPixel>
void VectorToImageFilter<TPixel>::Update() {
	m_tiler->Update();
}

} // End namespace itk

#endif // Define VECTORTOIMAGEFILTER_HXX
