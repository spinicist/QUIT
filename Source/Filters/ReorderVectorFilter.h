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
	void GenerateOutputInformation() ITK_OVERRIDE;

protected:
	ReorderVectorFilter() {}
	~ReorderVectorFilter() {}

    void ThreadedGenerateData(const TRegion &region, ThreadIdType threadId) ITK_OVERRIDE;
};

} // End namespace itk

#include "ReorderVectorFilter.hxx"

#endif // REORDERVECTORFILTER_H
