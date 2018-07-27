/*
 *  MeanImageFilter.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef MEANIMAGEFILTER_H
#define MEANIMAGEFILTER_H

#include "itkImageToImageFilter.h"

namespace itk {

template<typename TInput>
class MeanImageFilter : public ImageToImageFilter<TInput, TInput>
{
public:
    static const size_t ImageDimension = TInput::ImageDimension;
    typedef typename TInput::PixelType           TPixel;
    typedef typename TInput::RegionType          TRegion;
    typedef TInput                               TOutput;

    typedef MeanImageFilter                     Self;
    typedef ImageToImageFilter<TInput, TOutput> Superclass;
    typedef itk::SmartPointer<Self>             Pointer;

    itkNewMacro(Self);
    itkTypeMacro(MeanImageFilter, ImageToImageFilter);

    void SetNumberOfImages(size_t n);

protected:
    MeanImageFilter();
    ~MeanImageFilter(){}

    void GenerateOutputInformation() ITK_OVERRIDE; // Because output will be different dimension to input
    void ThreadedGenerateData(const TRegion &region, ThreadIdType threadId) ITK_OVERRIDE; // Does the work
    //DataObject::Pointer MakeOutput(unsigned int idx); // Create the Output

private:
    MeanImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

} // End namespace itk

#include "MeanImageFilter.hxx"

#endif // MEANIMAGEFILTER_H
