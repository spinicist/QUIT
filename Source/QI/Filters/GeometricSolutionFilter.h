/*
 *  GeometricSolutionFilter.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_GEOMETRICSOLUTIONFILTER_H
#define QI_GEOMETRICSOLUTIONFILTER_H

#include <complex>
#include <Eigen/Core>
#include "itkImageToImageFilter.h"

namespace QI {

// Helper function - complex equivalent of the dot product
template<typename T> inline T cdot(const std::complex<T> &a, const std::complex<T> &b) {
    return real(a * conj(b));
}

class GeometricSolutionFilter : public itk::ImageToImageFilter<itk::VectorImage<std::complex<float>, 3>, itk::Image<std::complex<float>, 3>> {
public:
    enum class RegEnum { None = 0, Line, Magnitude };
protected:
    size_t m_flips, m_lines, m_crossings, m_phases = 0;
    RegEnum m_Regularise = RegEnum::Line;
public:
    /** Standard class typedefs. */
    typedef itk::VectorImage<std::complex<float>, 3> TIn;
    typedef itk::Image<std::complex<float>, 3>       TOut;
    typedef itk::Image<float, 3>                TMask;
    typedef GeometricSolutionFilter             Self;
    typedef itk::ImageToImageFilter<TIn, TOut>  Superclass;
    typedef itk::SmartPointer<Self>             Pointer;

    itkNewMacro(Self); /** Method for creation through the object factory. */
    itkTypeMacro(GSFilter, ImageToImageFilter); /** Run-time type information (and related methods). */

    void SetRegularise(const RegEnum &r);
    const RegEnum &GetRegularise();
    void SetPhases(const size_t p);
    void SetMask(const TMask *mask);
    typename TMask::ConstPointer GetMask() const;
    virtual void GenerateOutputInformation() ITK_OVERRIDE;
protected:
    GeometricSolutionFilter();

    virtual void ThreadedGenerateData(const TIn::RegionType &region, itk::ThreadIdType threadId) ITK_OVERRIDE;
private:
    GeometricSolutionFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

} // End namespace QI

#include "QI/Filters/GeometricSolutionFilter.hxx"

#endif