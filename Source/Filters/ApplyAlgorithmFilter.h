#ifndef APPLYALGOFILTER_H
#define APPLYALGOFILTER_H

#include <vector>
#include "itkImageToImageFilter.h"
#include "itkVariableLengthVector.h"
#include "itkVectorImage.h"
#include "itkTimeProbe.h"
#include "QI/ThreadPool.h"

namespace itk{

template<typename TInputImage_, typename TOutputImage_, typename TConstImage_>
class ApplyAlgorithmFilter : public ImageToImageFilter<TInputImage_, TOutputImage_> {
public:
    typedef TInputImage_  TInputImage;
    typedef TOutputImage_ TOutputImage;
    typedef TConstImage_  TConstImage;

    typedef typename TInputImage::PixelType  TInputPixel;
    typedef typename TOutputImage::PixelType TOutputPixel;
    typedef typename TConstImage::PixelType  TConstPixel;
    typedef unsigned int TIterations;
    typedef Image<TIterations, TInputImage::ImageDimension> TIterationsImage;

    typedef ApplyAlgorithmFilter                          Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self>                            Pointer;
    typedef typename TInputImage::RegionType              TRegion;

    itkNewMacro(Self); /** Method for creation through the object factory. */
    itkTypeMacro(ApplyAlgorithmFilter, ImageToImageFilter); /** Run-time type information (and related methods). */

    class Algorithm {
    public:
        using TInput = TInputPixel;
        using TOutput = TOutputPixel;
        using TConst = TConstPixel;
        typedef TIterations  TIters;
        virtual size_t numInputs() const = 0;  // The number of inputs that will be concatenated into the data vector
        virtual size_t numConsts() const = 0;  // Number of constant input parameters/variables
        virtual size_t numOutputs() const = 0; // Number of output parameters/variables
        virtual size_t outputSize(const int i) const { return 1; }; // Size of each output, overload for vector outputs
        virtual size_t dataSize() const = 0;   // The expected size of the concatenated input
        virtual std::vector<TConst> defaultConsts() const = 0;    // Give some default constants for when the user does not supply them
        virtual void apply(const std::vector<TInput> &inputs,
                           const std::vector<TConst> &consts,
                           std::vector<TOutput> &outputs,
                           TConst &residual, TInput &resids,
                           TIters &iterations) const = 0; // Apply the algorithm to the data from one voxel
        virtual const TOutput &zero(const size_t i) const = 0; // Hack, to supply a zero for masked voxels
    };

	void SetAlgorithm(const std::shared_ptr<Algorithm> &a);
	std::shared_ptr<const Algorithm> GetAlgorithm() const;

	void SetInput(const size_t i, const TInputImage *img);
	typename TInputImage::ConstPointer GetInput(const size_t i) const;
	void SetConst(const size_t i, const TConstImage *img);
	typename TConstImage::ConstPointer GetConst(const size_t i) const;
	void SetMask(const TConstImage *mask);
	typename TConstImage::ConstPointer GetMask() const;

    void SetPoolsize(const size_t nThreads);
    void SetSubregion(const TRegion &sr); 
    void SetVerbose(const bool v);
    void SetOutputAllResiduals(const bool r); 
    
    TOutputImage     *GetOutput(const size_t i);
    TConstImage      *GetResidualOutput();
    TInputImage      *GetAllResidualsOutput();
    TIterationsImage *GetIterationsOutput();

    RealTimeClock::TimeStampType GetTotalTime() const;
    RealTimeClock::TimeStampType GetMeanTime() const;
    SizeValueType GetEvaluations() const;

protected:
	ApplyAlgorithmFilter();
	~ApplyAlgorithmFilter(){}
	DataObject::Pointer MakeOutput(unsigned int idx) ITK_OVERRIDE;

	std::shared_ptr<Algorithm> m_algorithm;
    bool m_verbose = false, m_hasSubregion = false, m_allResiduals = false;
    size_t m_poolsize = 1;
    TRegion m_subregion;

    RealTimeClock::TimeStampType m_elapsedTime = 0.0;
    SizeValueType m_unmaskedVoxels = 0;
    static const int AllResidualsOutputOffset = 0;
	static const int ResidualOutputOffset = 1;
    static const int IterationsOutputOffset = 2;
    static const int ExtraOutputs = 3;

    virtual void GenerateOutputInformation() ITK_OVERRIDE;
    virtual void GenerateData() ITK_OVERRIDE;

private:
	ApplyAlgorithmFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented

};

}

#include "ApplyAlgorithmFilter.hxx"

#endif // APPLYAGLOFILTER_H
