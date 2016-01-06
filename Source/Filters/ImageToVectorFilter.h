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
public:
	static const size_t InputDimension = TInput::ImageDimension;
	static const size_t OutputDimension = TInput::ImageDimension - 1;
	typedef typename TInput::PixelType           TPixel;
	typedef VectorImage<TPixel, OutputDimension> TOutput;
	typedef Image<TPixel, OutputDimension>       TVolume;

	typedef ImageToVectorFilter                 Self;
	typedef ImageToImageFilter<TInput, TOutput> Superclass;
	typedef itk::SmartPointer<Self>             Pointer;

	itkNewMacro(Self);
	itkTypeMacro(ImageToVectorFilter, ImageToImageFilter);

    itkSetMacro(BlockStart, size_t);
    itkSetMacro(BlockSize, size_t);

protected:
	typedef ExtractImageFilter<TInput, TVolume>  ExtractType;
	typedef ComposeImageFilter<TVolume, TOutput> ComposeType;

    typename ComposeType::Pointer m_compose;
    size_t m_BlockStart, m_BlockSize;

	ImageToVectorFilter();
	~ImageToVectorFilter(){}

    virtual void GenerateOutputInformation() override; // Because output will be different dimension to input
    virtual void GenerateData() override; // Does the work
	//DataObject::Pointer MakeOutput(unsigned int idx); // Create the Output

private:
	ImageToVectorFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
};

} // End namespace itk

#include "ImageToVectorFilter.hxx"

#endif // IMAGETOVECTORFILTER_H
