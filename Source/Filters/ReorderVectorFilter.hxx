#ifndef REORDERVECTORFILTER_HXX
#define REORDERVECTORFILTER_HXX

namespace itk {

template<typename TImage>
void ReorderVectorFilter<TImage>::GenerateOutputInformation() {
	//std::cout << __PRETTY_FUNCTION__ << std::endl;
	Superclass::GenerateOutputInformation();
	m_fullsize = this->GetInput()->GetNumberOfComponentsPerPixel();

	if (m_blocksize == 0)
		m_blocksize = m_fullsize;

	if ((m_fullsize % m_blocksize) != 0) {
		throw(std::runtime_error("Fullsize must be an integer multiple of blocksize."));
	}
	if ((m_blocksize % m_stride) != 0) {
		throw(std::runtime_error("Blocksize must be an integer multiple of stride."));
	}
	m_blocks = m_fullsize / m_blocksize;
	//std::cout << "End " << __PRETTY_FUNCTION__ << std::endl;
}

template<typename TImage>
void ReorderVectorFilter<TImage>::ThreadedGenerateData(const TRegion &region, ThreadIdType threadId) {
	//std::cout <<  __PRETTY_FUNCTION__ << std::endl;
	ImageRegionConstIterator<TImage> inIt(this->GetInput(), region);
	ImageRegionIterator<TImage>      outIt(this->GetOutput(), region);

	inIt.GoToBegin();
	outIt.GoToBegin();
	while(!inIt.IsAtEnd()) {
		VariableLengthVector<TPixel> in = inIt.Get();
		VariableLengthVector<TPixel> out(in.GetNumberOfElements());
		size_t outIndex = 0;
		for (size_t b = 0; b < m_fullsize; b += m_blocksize) { // Tracks the start of each block
			// Within each block, treat the data as matrix that needs transposing
			for (size_t j = b; j < (b + m_stride); j++) { // Tracks the 'column'
				for (size_t i = j; i < (b + m_blocksize); i += m_stride) { // Tracks the elements in 'rows'
					out[outIndex] = in[i];
					outIndex++;
				}
			}
		}
		outIt.Set(out);
		++inIt;
		++outIt;
	}
	// std::cout << "End " << __PRETTY_FUNCTION__ << std::endl;
}

} // End namespace itk

#endif // Define REORDERVECTORFILTER_HXX
