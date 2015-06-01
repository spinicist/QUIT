#ifndef IMAGETOVECTORFILTER_H
#define IMAGETOVECTORFILTER_H

#include "itkVectorImage.h"
#include "itkImageToImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkComposeImageFilter.h"

namespace itk {

template<typename TInput>
class ImageToVectorFilter : public ImageToImageFilter<TInput, VectorImage<typename TInput::PixelType, TInput::ImageDimension - 1>>
{
protected:
	size_t m_start, m_stop, m_stride = 0;

public:
	typedef typename TInput::PixelType                      TPixel;
	typedef VectorImage<TPixel, TInput::ImageDimension - 1> TOutput;
	typedef Image<TPixel, TInput::ImageDimension - 1>       TVolume;

	typedef ImageToVectorFilter                 Self;
	typedef ImageToImageFilter<TInput, TOutput> Superclass;
	typedef itk::SmartPointer<Self>             Pointer;

	itkNewMacro(Self);
	itkTypeMacro(ImageToVectorFilter, ImageToImageFilter);

	void SetStartStop(size_t start, size_t stop) const;
	void SetStride(size_t stride) const;

protected:
	typedef ExtractImageFilter<TInput, TVolume>  ExtractType;
	typedef ComposeImageFilter<TVolume, TOutput> ComposeType;
	typename ComposeType::Pointer m_compose;

	ImageToVectorFilter();
	~ImageToVectorFilter(){}

	virtual void GenerateOutputInformation(); // Because output will be different dimension to input
	virtual void GenerateData(); // Does the work
	//DataObject::Pointer MakeOutput(unsigned int idx); // Create the Output

private:
	ImageToVectorFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
};

} // End namespace itk

#include "ImageToVectorFilter.hxx"

#endif // IMAGETOVECTORFILTER_H
