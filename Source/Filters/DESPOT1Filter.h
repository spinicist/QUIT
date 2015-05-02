#ifndef DESPOT1FILTER_H
#define DESPOT1FILTER_H

#include "../Sequence.h"
#include "../Model.h"
#include "../Algorithm.h"

namespace itk{

template<typename TVectorImage, typename TImage>
class DESPOT1Filter : public ImageToImageFilter<TVectorImage, TImage>
{
public:
	/** Standard class typedefs. */
	typedef DESPOT1Filter                            Self;
	typedef ImageToImageFilter<TVectorImage, TImage> Superclass;
	typedef SmartPointer<Self>                       Pointer;
	typedef typename TImage::RegionType              RegionType;

	itkNewMacro(Self); /** Method for creation through the object factory. */
	itkTypeMacro(DESPOT1Filter, ImageToImageFilter); /** Run-time type information (and related methods). */

	void SetInput(const TVectorImage *SPGR);
	typename TVectorImage::ConstPointer GetInput();
	void SetMask(const TImage *mask);
	typename TImage::ConstPointer GetMask();
	void SetB1(const TImage *B1);
	typename TImage::ConstPointer GetB1();
	TImage *GetOutput(const size_t i);

	void SetSequence(const shared_ptr<SPGRSimple> &seq);
	void SetAlgorithm(const shared_ptr<Algorithm> &a);

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
	DESPOT1Filter();
	~DESPOT1Filter(){}

	virtual void ThreadedGenerateData(const RegionType & outputRegionForThread,
	                                  ThreadIdType threadId); // Does the work
	DataObject::Pointer MakeOutput(unsigned int idx); // Create the Output

private:
	DESPOT1Filter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented

	shared_ptr<SPGRSimple> m_sequence;
	shared_ptr<Algorithm> m_algorithm;
	size_t m_iterations;
};
}

#include "DESPOT1Filter.hxx"

#endif // DESPOT1FILTER_H
