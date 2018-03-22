/*
 *  VectorImageWrite.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <string>

#include "itkImageFileWriter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkDivideImageFilter.h"

#include "VectorToImageFilter.h"

#include "ImageIO.h"
#include "Macro.h"

namespace QI {

template<typename TVImg>
void WriteVectorImage(const TVImg *img, const std::string &path) {
    using TToSeries = itk::VectorToImageFilter<TVImg>;

    typename TToSeries::Pointer convert = TToSeries::New();
    convert->SetInput(img);
    convert->Update();
    WriteImage(convert->GetOutput(), path);
}

template<typename TVImg>
void WriteVectorImage(const itk::SmartPointer<TVImg> &ptr, const std::string &path) {
    WriteVectorImage(ptr.GetPointer(), path);
}

template<typename TVImg>
void WriteVectorMagnitudeImage(const TVImg *img, const std::string &path) {
    typedef itk::VectorToImageFilter<TVImg> TToSeries;
    auto convert = TToSeries::New();
    convert->SetInput(img);

    typedef typename TVImg::InternalPixelType TPixel;
    typedef typename TPixel::value_type TReal;
    typedef itk::Image<TPixel, 4> TSeries;
    typedef itk::Image<TReal, 4> TRealSeries;
    auto mag = itk::ComplexToModulusImageFilter<TSeries, TRealSeries>::New();
    mag->SetInput(convert->GetOutput());
    mag->Update();
    WriteImage<TRealSeries>(mag->GetOutput(), path);
}

template<typename TVImg>
void WriteVectorMagnitudeImage(const itk::SmartPointer<TVImg> &ptr, const std::string &path) {
    WriteVectorMagnitudeImage(ptr.GetPointer(), path);
}

template<typename TVImg>
void WriteScaledVectorImage(const TVImg *img, const QI::VolumeF *simg, const std::string &path) {
    auto scaleFilter = itk::DivideImageFilter<TVImg, QI::VolumeF, TVImg>::New();
    scaleFilter->SetInput1(img);
    scaleFilter->SetInput2(simg);
    scaleFilter->Update();
    WriteVectorImage(scaleFilter->GetOutput(), path);
}

template<typename TVImg>
void WriteScaledVectorImage(const itk::SmartPointer<TVImg> &ptr, const itk::SmartPointer<QI::VolumeF> &sptr, const std::string &path) {
    WriteScaledVectorImage(ptr.GetPointer(), sptr.GetPointer(), path);
}

template void WriteVectorImage<VectorVolumeF>(const VectorVolumeF *img, const std::string &path);
template void WriteVectorImage<VectorVolumeXF>(const VectorVolumeXF *img, const std::string &path);
template void WriteVectorImage<VectorVolumeF>(const itk::SmartPointer<VectorVolumeF> &ptr, const std::string &path);
template void WriteVectorImage<VectorVolumeXF>(const itk::SmartPointer<VectorVolumeXF> &ptr, const std::string &path);
template void WriteVectorMagnitudeImage<VectorVolumeXF>(const VectorVolumeXF *ptr, const std::string &path);
template void WriteVectorMagnitudeImage<VectorVolumeXF>(const itk::SmartPointer<VectorVolumeXF> &ptr, const std::string &path);
template void WriteScaledVectorImage<VectorVolumeF>(const VectorVolumeF *img, const VolumeF *simg, const std::string &path);
template void WriteScaledVectorImage<VectorVolumeF>(const itk::SmartPointer<VectorVolumeF> &ptr, const itk::SmartPointer<VolumeF> &sptr, const std::string &path);
} // End namespace QUIT
