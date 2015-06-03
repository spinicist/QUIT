#ifndef VECTORTOIMAGEFILTER_H
#define VECTORTOIMAGEFILTER_H

#include "itkImageToImageFilter.h"
#include "itkVectorImage.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkTileImageFilter.h"

namespace itk {

template<typename TPixel>
class VectorToImageFilter  : public itk::ImageToImageFilter<itk::VectorImage<TPixel, 3>, itk::Image<TPixel, 4>>
{
public:
	typedef itk::Image<TPixel, 3> TImage;
	typedef itk::VectorImage<TPixel, 3> TInput;
	typedef itk::Image<TPixel, 4> TOutput;

	typedef VectorToImageFilter                            Self;
	typedef itk::ImageToImageFilter<TInput, TOutput> Superclass;
	typedef itk::SmartPointer<Self>                        Pointer;

	itkNewMacro(Self);
	itkTypeMacro(Self, Superclass);

protected:
	typedef itk::VectorIndexSelectionCastImageFilter<TInput, TImage> TIndexer;
	typedef itk::TileImageFilter<TImage, TOutput> TTiler;
	typename TTiler::Pointer m_tiler;
	std::vector<typename TIndexer::Pointer> m_indexers;

	VectorToImageFilter();
	~VectorToImageFilter(){}

	virtual void GenerateOutputInformation(); // Because output will be different dimension to input
	virtual void GenerateData() {} // Don't need a GenerateData, just call m_tiler->Update()
	virtual void Update();
	//DataObject::Pointer MakeOutput(unsigned int idx); // Create the Output

private:
	VectorToImageFilter(const Self &); //purposely not implemented
	void operator=(const Self &);      //purposely not implemented
};

} // End namespace itk

#include "VectorToImageFilter.hxx"

#endif // VECTORTOIMAGEFILTER_H
