#ifndef VECTORTOIMAGEFILTER_HXX
#define VECTORTOIMAGEFILTER_HXX

namespace itk {

template<typename TInput>
VectorToImageFilter<TInput>::VectorToImageFilter() {
	m_tiler = TTiler::New();
}

template<typename TInput>
void VectorToImageFilter<TInput>::SetInput(const TInput *img) {
    Superclass::SetInput(img);
    
    auto in  = this->GetInput();
    if (!in) {
        return;
    }
    itk::FixedArray<unsigned int, 4> layout;
    layout[0] = layout[1] = layout[2] = 1;
    layout[3] = in->GetNumberOfComponentsPerPixel();
    m_tiler->SetLayout(layout);
    m_indexers.resize(in->GetNumberOfComponentsPerPixel());
    for (size_t i =0; i < in->GetNumberOfComponentsPerPixel(); i++) {
        m_indexers[i] = TIndexer::New();
        m_indexers[i]->SetIndex(i);
        m_indexers[i]->SetInput(in);
        m_tiler->SetInput(i, m_indexers[i]->GetOutput());
    }
    m_tiler->UpdateLargestPossibleRegion();
    this->GraftOutput(m_tiler->GetOutput());
}

template<typename TInput>
void VectorToImageFilter<TInput>::GenerateOutputInformation() {
    auto out = this->GetOutput();
    auto in  = this->GetInput();
    if ( !out || !in ) {
        return;
    }
    auto spacing    = in->GetSpacing();
    auto origin     = in->GetOrigin();
    auto direction  = in->GetDirection();
    auto inRegion   = in->GetLargestPossibleRegion();
    auto outSpacing   = out->GetSpacing();
    auto outOrigin    = out->GetOrigin();
    auto outDirection = out->GetDirection();
    typename TOutput::SizeType outSize;
    typename TOutput::IndexType outIndex;
    typename TOutput::RegionType outRegion;
    
    outSpacing.Fill(1);
    outOrigin.Fill(1);
    outDirection.SetIdentity();
    for (int i = 0; i < InputDimension; i++) {
        outSpacing[i] = spacing[i];
        outOrigin[i] =  origin[i];
        outSize[i] = inRegion.GetSize()[i];
        outIndex[i] = inRegion.GetIndex()[i];
        for (int j = 0; j < InputDimension; j++) {
            outDirection[i][j] = direction[i][j];
        }
    }
    out->SetSpacing(outSpacing);
    out->SetOrigin(outOrigin);
    out->SetDirection(outDirection);
    outIndex[OutputDimension - 1] = 0;
    outSize[OutputDimension - 1] = in->GetNumberOfComponentsPerPixel();
    outRegion.SetSize(outSize);
    outRegion.SetIndex(outIndex);
    out->SetLargestPossibleRegion(outRegion);
}

template<typename TInput>
void VectorToImageFilter<TInput>::GenerateData() {
    typename TInput::Pointer input = TInput::New();
    input->Graft(const_cast<TInput *>(this->GetInput()));
    auto spacing = this->GetOutput()->GetSpacing();
    auto direction = this->GetOutput()->GetDirection();
    auto origin = this->GetOutput()->GetOrigin();
    m_tiler->UpdateLargestPossibleRegion();
    this->GraftOutput(m_tiler->GetOutput());
    // Reset space information because tiler messes it up
    this->GetOutput()->SetSpacing(spacing);
    this->GetOutput()->SetDirection(direction);
    this->GetOutput()->SetOrigin(origin);
}

} // End namespace itk

#endif // Define VECTORTOIMAGEFILTER_HXX
