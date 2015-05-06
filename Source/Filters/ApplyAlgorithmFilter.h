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
		virtual size_t numConsts() const = 0;
		virtual size_t numOutputs() const = 0;

		virtual void apply(const shared_ptr<SequenceBase> sequence,
		                   const TInputVector &data,
		                   const VectorXd &consts,
		                   VectorXd &outputs,
		                   ArrayXd &resids) const = 0;

		virtual VectorXd defaultConsts() = 0;
};

namespace itk{

template<typename TData, typename TAlgo>
class ApplyAlgorithmFilter : public ImageToImageFilter<VectorImage<TData, 3>, Image<float, 3>>
{
public:
	/** Standard class typedefs. */
	typedef Image<float, 3>                         TImage;
	typedef VectorImage<TData, 3>                   TInputImage;
	typedef VectorImage<float, 3>                   TResidImage;
	typedef ApplyAlgorithmFilter                    Self;
	typedef ImageToImageFilter<TInputImage, TImage> Superclass;
	typedef SmartPointer<Self>                      Pointer;
	typedef typename TImage::RegionType             RegionType;

	itkNewMacro(Self); /** Method for creation through the object factory. */
	itkTypeMacro(ApplyAlgorithmFilter, ImageToImageFilter); /** Run-time type information (and related methods). */

	void SetDataInput(const size_t i, const TInputImage *img);
	void SetConstInput(const size_t i, const TImage *img);
	void SetMask(const TImage *mask);
	typename TInputImage::ConstPointer GetDataInput(const size_t i) const;
	typename TImage::ConstPointer GetConstInput(const size_t i) const;
	typename TImage::ConstPointer GetMask() const;
	TImage *GetOutput(const size_t i);
	TResidImage *GetResidOutput();

	void SetSequence(const shared_ptr<SequenceBase> &seq);
	void SetAlgorithm(const shared_ptr<TAlgo> &a);
	void Setup();

	virtual void GenerateOutputInformation() override;
	virtual void Update() override;

protected:
	ApplyAlgorithmFilter();
	~ApplyAlgorithmFilter(){}

	virtual void ThreadedGenerateData(const RegionType & outputRegionForThread,
	                                  ThreadIdType threadId); // Does the work
	DataObject::Pointer MakeOutput(unsigned int idx); // Create the Output

	shared_ptr<SequenceBase> m_sequence;
	shared_ptr<TAlgo> m_algorithm;

private:
	ApplyAlgorithmFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
};
}

#include "ApplyAlgorithmFilter.hxx"

#endif // APPLYAGLOFILTER_H
