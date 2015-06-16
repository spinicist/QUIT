#ifndef APPLYALGOSLICEBYSLICEFILTER_H
#define APPLYALGOSLICEBYSLICEFILTER_H

#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "ApplyAlgorithmFilter.h"

namespace itk{

template<typename TVImage, typename TAlgo>
class ApplyAlgorithmSliceBySliceFilter : public ApplyAlgorithmFilter<TVImage, TAlgo> {
public:
	typedef ApplyAlgorithmSliceBySliceFilter       Self;
	typedef ApplyAlgorithmFilter<TVImage, TVImage> Superclass;
	typedef SmartPointer<Self>                     Pointer;

	static const unsigned int                   ImageDimension = Superclass::ImageDimension;
	typedef typename Superclass::TPixel         TPixel;
	typedef typename Superclass::TImage         TImage;
	typedef typename Superclass::TVectorImage   TVectorImage;

	static const unsigned int                   SliceDimension = Superclass::ImageDimension - 1;
	typedef Image<TPixel, SliceDimension>       TSlice;
	typedef VectorImage<TPixel, SliceDimension> TVectorSlice;
	typedef ApplyAlgorithmFilter<TVectorSlice, TAlgo> TSliceFilter;

	itkNewMacro(Self);
	itkTypeMacro(ApplyAlgorithmSliceBySliceFilter, ApplyAlgorithmFilter);

	virtual void PrintSelf(std::ostream & os, Indent indent) const;

protected:
	ApplyAlgorithmSliceBySliceFilter();
	~ApplyAlgorithmSliceBySliceFilter(){}

	int m_sliceIndex = 0;

	void VerifyInputInformation() override;
	//void GenerateInputRequestedRegion() override;
	void GenerateData() override;
	DataObject::Pointer MakeOutput(unsigned int idx);

private:
	ApplyAlgorithmSliceBySliceFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
};

} // End namespace itk

#include "ApplyAlgorithmSliceBySliceFilter.hxx"

#endif // APPLYALGOSLICEBYSLICEFILTER_H
