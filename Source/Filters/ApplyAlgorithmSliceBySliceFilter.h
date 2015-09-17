#ifndef APPLYALGOSLICEBYSLICEFILTER_H
#define APPLYALGOSLICEBYSLICEFILTER_H

#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "ApplyAlgorithmFilter.h"

namespace itk{

template<typename TAlgorithm , typename TData = float, typename TScalar = float, unsigned int ImageDim = 3>
class ApplyAlgorithmSliceBySliceFilter :
	public ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim> {
public:
	typedef ApplyAlgorithmSliceBySliceFilter                           Self;
	typedef ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, ImageDim> Superclass;
	typedef SmartPointer<Self>                                         Pointer;

	typedef typename Superclass::TDataVectorImage   TDataVectorImage;
	typedef typename Superclass::TScalarVectorImage TScalarVectorImage;
	typedef typename Superclass::TScalarImage       TScalarImage;

	static const unsigned int SliceDim =   ImageDim - 1;
	typedef VectorImage<TData, SliceDim>   TDataVectorSlice;
	typedef VectorImage<TScalar, SliceDim> TScalarVectorSlice;
	typedef Image<TScalar, SliceDim>       TScalarSlice;

	typedef ApplyAlgorithmFilter<TAlgorithm, TData, TScalar, SliceDim> TSliceFilter;

	itkNewMacro(Self);
	itkTypeMacro(ApplyAlgorithmSliceBySliceFilter, ApplyAlgorithmFilter);

    virtual void PrintSelf(std::ostream & os, Indent indent) const override;
	void SetSlices(const int start, const int stop);
	int GetSliceIndex() const;

protected:
	ApplyAlgorithmSliceBySliceFilter();
	~ApplyAlgorithmSliceBySliceFilter(){}

	int m_sliceIndex = 0, m_startSlice = 0, m_stopSlice = 0;

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
