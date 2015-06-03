#ifndef REORDERVECTORFILTER_H
#define REORDERVECTORFILTER_H

#include "itkInPlaceImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkVectorImage.h"

namespace itk {

template<typename TImage>
class ReorderVectorFilter : public InPlaceImageFilter<TImage> {
protected:
	size_t m_stride = 1, m_fullsize = 0, m_blocksize = 0, m_blocks = 0;

public:
	typedef ReorderVectorFilter         Self;
	typedef InPlaceImageFilter<TImage>  Superclass;
	typedef SmartPointer<Self>          Pointer;
	typedef typename TImage::RegionType TRegion;
	typedef typename TImage::InternalPixelType TPixel;

	itkNewMacro(Self);
	itkTypeMacro(ReorderVectorFilter, InPlaceImageFilter);

	void SetStride(size_t s) { m_stride = s; }
	void SetBlockSize(size_t b) { m_blocksize = b; }

private:
	ReorderVectorFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented

public:
	virtual void GenerateOutputInformation() override {
		//std::cout << __PRETTY_FUNCTION__ << std::endl;
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
		Superclass::GenerateOutputInformation();
		//std::cout << "End " << __PRETTY_FUNCTION__ << std::endl;
	}

protected:
	ReorderVectorFilter() {}
	~ReorderVectorFilter() {}

	virtual void ThreadedGenerateData(const TRegion &region, ThreadIdType threadId) {
		//std::cout <<  __PRETTY_FUNCTION__ << std::endl;
		ImageRegionConstIterator<TImage> inIt(this->GetInput(), region);
		ImageRegionIterator<TImage>      outIt(this->GetOutput(), region);

		// A stride of 1 means do nothing
		if (m_stride == 1)
			return;

		while(!inIt.IsAtEnd()) {
			VariableLengthVector<TPixel> in = inIt.Get();
			VariableLengthVector<TPixel> out(in.GetNumberOfElements());

			size_t outIndex = 0;
			//std::cout << "blocksize " << m_blocksize << " blocks " << m_blocks << " fullsize " << m_fullsize << " stride " << m_stride << std::endl;
			for (size_t b = 0; b < m_fullsize; b += m_blocksize) { // Tracks the start of each block
				// Within each block, treat the data as matrix that needs transposing
				for (size_t j = b; j < m_stride; j++) { // Tracks the 'column'
					for (size_t i = j; i < m_blocksize; i += m_stride) { // Tracks the elements in 'rows'
						//std::cout << "b " << b << " j " << j << " i " << i << " out " << outIndex << std::endl;
						out[outIndex] = in[i];
						outIndex++;
					}
				}
			}

			//std::cout << "IN:  " << in << std::endl << "OUT: " << out << std::endl;
			outIt.Set(out);
			++inIt;
			++outIt;
		}
		//std::cout << "End " << __PRETTY_FUNCTION__ << std::endl;
	}
};

} // End namespace itk

#include "ReorderVectorFilter.hxx"

#endif // VECTORTOIMAGEFILTER_H
