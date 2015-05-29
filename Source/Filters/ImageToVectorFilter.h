#ifndef IMAGETOVECTORFILTER_H
#define IMAGETOVECTORFILTER_H

#include "itkVectorImage.h"
#include "itkImageToImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkComposeImageFilter.h"

template<typename PixelType>
class ImageToVectorFilter : public itk::ImageToImageFilter<itk::Image<PixelType, 4>, itk::VectorImage<PixelType, 3>>
{
public:
	typedef itk::Image<PixelType, 4> InputType;
	typedef itk::Image<PixelType, 3> VolumeType;
	typedef itk::VectorImage<PixelType, 3> OutputType;

	typedef ImageToVectorFilter                            Self;
	typedef itk::ImageToImageFilter<InputType, OutputType> Superclass;
	typedef itk::SmartPointer<Self>                        Pointer;

	itkNewMacro(Self);
	itkTypeMacro(ImageToVectorFilter, ImageToImageFilter);

protected:
	typedef itk::ExtractImageFilter<InputType, VolumeType>  ExtractType;
	typedef itk::ComposeImageFilter<VolumeType, OutputType> ComposeType;
	//typename ExtractType::Pointer m_ExtractVolume
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


#include "ImageToVectorFilter.hxx"

#endif // IMAGETOVECTORFILTER_H
