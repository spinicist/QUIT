#ifndef APPLYALGOFILTER_H
#define APPLYALGOFILTER_H

#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "../Sequence.h"
#include "../Model.h"

class Algorithm {
	public:
		virtual size_t numConsts() const = 0;
		virtual size_t numOutputs() const = 0;

		virtual void apply(const shared_ptr<SteadyState> sequence,
		                   const VectorXd &data,
		                   const VectorXd &consts,
		                   VectorXd &outputs,
		                   ArrayXd &resids) const = 0;

		virtual VectorXd defaultConsts() = 0;
};

namespace itk{

template<typename TVectorImage, typename TImage>
class ApplyAlgorithmFilter : public ImageToImageFilter<TVectorImage, TImage>
{
public:
	/** Standard class typedefs. */
	typedef ApplyAlgorithmFilter                            Self;
	typedef ImageToImageFilter<TVectorImage, TImage> Superclass;
	typedef SmartPointer<Self>                       Pointer;
	typedef typename TImage::RegionType              RegionType;

	itkNewMacro(Self); /** Method for creation through the object factory. */
	itkTypeMacro(ApplyAlgorithmFilter, ImageToImageFilter); /** Run-time type information (and related methods). */

	void SetDataInput(const size_t i, const TVectorImage *img);
	void SetConstInput(const size_t i, const TImage *img);
	void SetMask(const TImage *mask);
	typename TVectorImage::ConstPointer GetDataInput(const size_t i) const;
	typename TImage::ConstPointer GetConstInput(const size_t i) const;
	typename TImage::ConstPointer GetMask() const;
	TImage *GetOutput(const size_t i);
	TVectorImage *GetResidOutput();

	void SetSequence(const shared_ptr<SPGRSimple> &seq);
	void SetAlgorithm(const shared_ptr<Algorithm> &a);
	void Setup();

	virtual void GenerateOutputInformation();
	virtual void Update();

	void PrintDirections() {
		cout << __PRETTY_FUNCTION__ << endl;
		typedef ImageBase< 3 > ImageBaseType;
		ImageBaseType *ptr = ITK_NULLPTR;
		InputDataObjectIterator it(this);

		for(; !it.IsAtEnd(); ++it ) {
			ptr = dynamic_cast< ImageBaseType * >( it.GetInput() );
			cout << it.GetName() << endl << ptr->GetDirection() << endl;
		}
	}

protected:
	ApplyAlgorithmFilter();
	~ApplyAlgorithmFilter(){}

	virtual void ThreadedGenerateData(const RegionType & outputRegionForThread,
	                                  ThreadIdType threadId); // Does the work
	DataObject::Pointer MakeOutput(unsigned int idx); // Create the Output

	shared_ptr<SPGRSimple> m_sequence;
	shared_ptr<Algorithm> m_algorithm;

private:
	ApplyAlgorithmFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
};
}

#include "ApplyAlgorithmFilter.hxx"

#endif // APPLYAGLOFILTER_H
