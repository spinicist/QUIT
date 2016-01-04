#ifndef REORDERIMAGEFILTER_H
#define REORDERIMAGEFILTER_H

#include "itkInPlaceImageFilter.h"

namespace itk {

template<typename TImage>
class ReorderImageFilter : public ImageToImageFilter<TImage, TImage> {
protected:
    size_t m_stride = 1, m_fullsize = 0, m_blocksize = 0, m_blocks = 0;

public:
    typedef ReorderImageFilter          Self;
    typedef InPlaceImageFilter<TImage>  Superclass;
    typedef SmartPointer<Self>          Pointer;
    typedef typename TImage::RegionType TRegion;
    typedef typename TImage::PixelType  TPixel;

    itkNewMacro(Self);
    itkTypeMacro(ReorderImageFilter, InPlaceImageFilter);

    void SetStride(size_t s) { m_stride = s; }
    void SetBlockSize(size_t b) { m_blocksize = b; }

private:
    ReorderImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented

public:
    virtual void GenerateOutputInformation() override;

protected:
    ReorderImageFilter() {}
    ~ReorderImageFilter() {}

    virtual void GenerateData() override;

};

} // End namespace itk

#include "ReorderImageFilter.hxx"

#endif // REORDERIMAGEFILTER_H
