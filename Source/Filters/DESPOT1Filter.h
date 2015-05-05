#ifndef DESPOT1FILTER_H
#define DESPOT1FILTER_H

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
		                   const VectorXd &inputs,
		                   VectorXd &outputs,
		                   ArrayXd &resids) const = 0;
};

class D1 : public Algorithm {
	public:
		enum class Type { LLS, WLLS, NLLS };
	private:
		Type m_type = Type::LLS;
		size_t m_iterations = 4;

		// T1 only Functor
class T1Functor : public DenseFunctor<double> {
	protected:
		const shared_ptr<SteadyState> m_sequence;
		const ArrayXd m_data;
		const bool m_debug;
		const double m_B1;
		const shared_ptr<SCD> m_model = make_shared<SCD>();

	public:
		T1Functor(const shared_ptr<SteadyState> cs, const ArrayXd &data,
		          const double B1, const bool debug) :
			DenseFunctor<double>(2, cs->size()),
			m_sequence(cs), m_data(data),
			m_B1(B1), m_debug(debug)
		{
			assert(static_cast<size_t>(m_data.rows()) == values());
		}

		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
			eigen_assert(diffs.size() == values());
			VectorXd fullParams = VectorXd::Zero(5);
			fullParams.head(2) = params;
			fullParams(4) = m_B1;
			ArrayXcd s = m_sequence->signal(m_model, fullParams);
			diffs = s.abs() - m_data;
			if (m_debug) {
				cout << endl << __PRETTY_FUNCTION__ << endl;
				cout << "p:     " << params.transpose() << endl;
				cout << "s:     " << s.transpose() << endl;
				cout << "data:  " << m_data.transpose() << endl;
				cout << "diffs: " << diffs.transpose() << endl;
			}
			return 0;
		}
};

	public:
		void setType(Type t) { m_type = t; }
		void setIterations(size_t n) { m_iterations = n; }

		size_t numConsts() const override { return 1; }
		size_t numOutputs() const override { return 2; }

		virtual void apply(const shared_ptr<SteadyState> sequence,
		                   const VectorXd &data,
		                   const VectorXd &inputs,
		                   VectorXd &outputs,
		                   ArrayXd &resids) const override
		{
			double B1 = inputs[0];
			ArrayXd flip = sequence->flip() * B1;
			VectorXd Y = data.array() / flip.sin();
			MatrixXd X(Y.rows(), 2);
			X.col(0) = data.array() / flip.tan();
			X.col(1).setOnes();
			VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
			outputs[1] = -sequence->TR() / log(b[0]);
			outputs[0] = b[1] / (1. - b[0]);
			if (m_type == Type::WLLS) {
				VectorXd W(sequence->size());
				for (size_t n = 0; n < m_iterations; n++) {
					W = (flip.sin() / (1. - (exp(-sequence->TR()/outputs[1])*flip.cos()))).square();
					b = (X.transpose() * W.asDiagonal() * X).partialPivLu().solve(X.transpose() * W.asDiagonal() * Y);
					outputs[1] = -sequence->TR() / log(b[0]);
					outputs[0] = b[1] / (1. - b[0]);
				}
			} else if (m_type == Type::NLLS) {
				T1Functor f(sequence, data, B1, false);
				NumericalDiff<T1Functor> nDiff(f);
				LevenbergMarquardt<NumericalDiff<T1Functor>> lm(nDiff);
				lm.setMaxfev(m_iterations * (sequence->size() + 1));
				lm.minimize(outputs);
			}

			ArrayXd theory = One_SPGR(sequence->flip(), sequence->TR(), outputs[0], outputs[1], B1).array().abs();
			resids = (data.array() - theory);
		}
};

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
