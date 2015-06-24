#ifndef APPLYALGOFILTER_H
#define APPLYALGOFILTER_H

#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "../Sequence.h"
#include "../Model.h"

#include "itkImageToImageFilter.h"
#include "itkSliceBySliceImageFilter.h"

template<typename DataType>
class Algorithm {
public:
	typedef DataType TScalar;
	typedef Eigen::Array<TScalar, Eigen::Dynamic, 1> TInput;
	typedef Eigen::ArrayXd TArray;
	virtual size_t numInputs() const = 0;  // The number of inputs that will be concatenated into the data vector
	virtual size_t numConsts() const = 0;  // Number of constant input parameters/variables
	virtual size_t numOutputs() const = 0; // Number of output parameters/variables
	virtual size_t dataSize() const = 0;   // The expected size of the concatenated data vector

	virtual void apply(const TInput &data,
					   const TArray &consts,
					   TArray &outputs,
					   TArray &resids) const = 0;

	virtual TArray defaultConsts() = 0;
};

namespace itk{

template<typename TInImage, typename TAlgo>
class ApplyAlgorithmFilter : public ImageToImageFilter<TInImage, TInImage> {
public:
	static const unsigned int                    ImageDimension = TInImage::ImageDimension;
	typedef typename TInImage::InternalPixelType TPixel;
	typedef typename TInImage::PixelType         TVector;
	typedef Image<TPixel, ImageDimension>        TImage;
	typedef TInImage                             TVectorImage;
	typedef TAlgo                                TAlgorithm;

	typedef ApplyAlgorithmFilter                   Self;
	typedef ImageToImageFilter<TInImage, TInImage> Superclass;
	typedef SmartPointer<Self>                     Pointer;
	typedef typename TImage::RegionType            TRegion;

	itkNewMacro(Self); /** Method for creation through the object factory. */
	itkTypeMacro(ApplyAlgorithmFilter, ImageToImageFilter); /** Run-time type information (and related methods). */

	void SetAlgorithm(const shared_ptr<TAlgo> &a);
	shared_ptr<const TAlgo> GetAlgorithm() const;
	void SetScaleToMean(const bool s);
	bool GetScaleToMean() const;

	void SetInput(const size_t i, const TVectorImage *img);
	typename TVectorImage::ConstPointer GetInput(const size_t i) const;
	void SetConst(const size_t i, const TImage *img);
	typename TImage::ConstPointer GetConst(const size_t i) const;
	void SetMask(const TImage *mask);
	typename TImage::ConstPointer GetMask() const;

	TImage *GetOutput(const size_t i);
	TVectorImage *GetResidOutput();

	virtual void GenerateOutputInformation() override;

protected:
	ApplyAlgorithmFilter();
	~ApplyAlgorithmFilter(){}

	virtual void ThreadedGenerateData(const TRegion &outputRegionForThread, ThreadIdType threadId) override;
	DataObject::Pointer MakeOutput(unsigned int idx);

	shared_ptr<TAlgo> m_algorithm;
	bool m_scale_to_mean = false;

private:
	ApplyAlgorithmFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
};
}

#include "ApplyAlgorithmFilter.hxx"

#endif // APPLYAGLOFILTER_H
