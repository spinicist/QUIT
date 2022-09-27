#ifndef VECTORTOIMAGEFILTER_H
#define VECTORTOIMAGEFILTER_H

#include "itkImageToImageFilter.h"
#include "itkVectorImage.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkTileImageFilter.h"

namespace itk {

template<typename TInput>
class VectorToImageFilter  : public ImageToImageFilter<TInput, Image<typename TInput::InternalPixelType, TInput::ImageDimension + 1>>
{
public:
	static const size_t InputDimension = TInput::ImageDimension;
	static const size_t OutputDimension = TInput::ImageDimension + 1;
	typedef typename TInput::InternalPixelType TPixel;
	typedef Image<TPixel, OutputDimension>     TOutput;
	typedef Image<TPixel, InputDimension>      TVolume;


	typedef VectorToImageFilter                      Self;
	typedef itk::ImageToImageFilter<TInput, TOutput> Superclass;
	typedef itk::SmartPointer<Self>                  Pointer;

	itkNewMacro(Self)
	itkTypeMacro(Self, Superclass);

    void SetInput(const TInput *img) ITK_OVERRIDE;
    
protected:
	typedef itk::VectorIndexSelectionCastImageFilter<TInput, TVolume> TIndexer;
	typedef itk::TileImageFilter<TVolume, TOutput> TTiler;
	typename TTiler::Pointer m_tiler;
	std::vector<typename TIndexer::Pointer> m_indexers;

	VectorToImageFilter();
	~VectorToImageFilter(){}

    void GenerateOutputInformation() ITK_OVERRIDE; // Because output will be different dimension to input
    void GenerateData() ITK_OVERRIDE; // Don't need a GenerateData, just call m_tiler->Update()

private:
	VectorToImageFilter(const Self &); //purposely not implemented
	void operator=(const Self &);      //purposely not implemented
};

} // End namespace itk

#include "VectorToImageFilter.hxx"

#endif // VECTORTOIMAGEFILTER_H
