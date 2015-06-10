#ifndef APPLYALGOFILTER_H
#define APPLYALGOFILTER_H

#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "../Sequence.h"
#include "../Model.h"

template<typename DataType>
class Algorithm {
public:
	typedef DataType TInput;
	typedef Eigen::Matrix<DataType, Eigen::Dynamic, 1> TInputVector;
	typedef Eigen::VectorXd TConstVector;
	virtual size_t numInputs() const = 0;  // The number of inputs that will be concatenated into the data vector
	virtual size_t numConsts() const = 0;  // Number of constant input parameters/variables
	virtual size_t numOutputs() const = 0; // Number of output parameters/variables
	virtual size_t dataSize() const = 0;   // The expected size of the concatenated data vector

	virtual void apply(const TInputVector &data,
					   const TConstVector &consts,
					   VectorXd &outputs,
					   ArrayXd &resids) const = 0;

	virtual TConstVector defaultConsts() = 0;
};

namespace itk{

template<typename TInImage, typename TAlgo>
class ApplyAlgorithmFilter : public ImageToImageFilter<TInImage, TInImage> {
public:
	static const unsigned int                    ImageDimension = TInImage::ImageDimension;
	typedef typename TInImage::InternalPixelType TPixel;
	typedef Image<TPixel, ImageDimension>        TImage;
	typedef TInImage                             TVectorImage;

	typedef ApplyAlgorithmFilter                   Self;
	typedef ImageToImageFilter<TInImage, TInImage> Superclass;
	typedef SmartPointer<Self>                     Pointer;
	typedef typename TImage::RegionType            TRegion;

	itkNewMacro(Self); /** Method for creation through the object factory. */
	itkTypeMacro(ApplyAlgorithmFilter, ImageToImageFilter); /** Run-time type information (and related methods). */

	void SetAlgorithm(const shared_ptr<TAlgo> &a);
	void SetDataInput(const size_t i, const TVectorImage *img);
	void SetConstInput(const size_t i, const TImage *img);
	void SetMask(const TImage *mask);
	void SetSlices(const int start, const int stop);
	typename TVectorImage::ConstPointer GetDataInput(const size_t i) const;
	typename TImage::ConstPointer GetConstInput(const size_t i) const;
	typename TImage::ConstPointer GetMask() const;
	TImage *GetOutput(const size_t i);
	TVectorImage *GetResidOutput();

	virtual void GenerateOutputInformation() override;
	virtual void Update() override;

protected:
	ApplyAlgorithmFilter();
	~ApplyAlgorithmFilter(){}

	virtual void ThreadedGenerateData(const TRegion &outputRegionForThread, ThreadIdType threadId);
	DataObject::Pointer MakeOutput(unsigned int idx);

	shared_ptr<TAlgo> m_algorithm;
	int m_startSlice = 0, m_stopSlice = 0;

private:
	ApplyAlgorithmFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
};
}

#include "ApplyAlgorithmFilter.hxx"

#endif // APPLYAGLOFILTER_H
