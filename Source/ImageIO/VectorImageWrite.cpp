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

#include "itkComplexToModulusImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkImageFileWriter.h"

#include "VectorToImageFilter.h"

#include "ImageIO.h"
#include "Log.h"

namespace QI {

template <typename TVImg>
void WriteVectorImage(const TVImg *img, const std::string &path, const bool verbose) {
    using TToSeries = itk::VectorToImageFilter<TVImg>;

    typename TToSeries::Pointer convert = TToSeries::New();
    convert->SetInput(img);
    convert->Update();
    QI::Log(verbose, "Writing image: {}", path);
    WriteImage(convert->GetOutput(), path, verbose);
}

template <typename TVImg>
void WriteVectorImage(const itk::SmartPointer<TVImg> &ptr, const std::string &path,
                      const bool verbose) {
    WriteVectorImage(ptr.GetPointer(), path, verbose);
}

template <typename TVImg>
void WriteVectorMagnitudeImage(const TVImg *img, const std::string &path, const bool verbose) {
    typedef itk::VectorToImageFilter<TVImg> TToSeries;
    auto                                    convert = TToSeries::New();
    convert->SetInput(img);

    typedef typename TVImg::InternalPixelType TPixel;
    typedef typename TPixel::value_type       TReal;
    typedef itk::Image<TPixel, 4>             TSeries;
    typedef itk::Image<TReal, 4>              TRealSeries;
    auto mag = itk::ComplexToModulusImageFilter<TSeries, TRealSeries>::New();
    mag->SetInput(convert->GetOutput());
    mag->Update();
    WriteImage<TRealSeries>(mag->GetOutput(), path, verbose);
}

template <typename TVImg>
void WriteVectorMagnitudeImage(const itk::SmartPointer<TVImg> &ptr, const std::string &path,
                               const bool verbose) {
    WriteVectorMagnitudeImage(ptr.GetPointer(), path, verbose);
}

template <typename TVImg>
void WriteScaledVectorImage(const TVImg *img, const QI::VolumeF *simg, const std::string &path,
                            const bool verbose) {
    auto scaleFilter = itk::DivideImageFilter<TVImg, QI::VolumeF, TVImg>::New();
    scaleFilter->SetInput1(img);
    scaleFilter->SetInput2(simg);
    scaleFilter->Update();
    WriteVectorImage(scaleFilter->GetOutput(), path, verbose);
}

template <typename TVImg>
void WriteScaledVectorImage(const itk::SmartPointer<TVImg> &      ptr,
                            const itk::SmartPointer<QI::VolumeF> &sptr, const std::string &path,
                            const bool verbose) {
    WriteScaledVectorImage(ptr.GetPointer(), sptr.GetPointer(), path, verbose);
}

template void WriteVectorImage<VectorVolumeF>(const VectorVolumeF *img, const std::string &path,
                                              const bool verbose);
template void WriteVectorImage<VectorVolumeXF>(const VectorVolumeXF *img, const std::string &path,
                                               const bool verbose);
template void WriteVectorImage<VectorVolumeF>(const itk::SmartPointer<VectorVolumeF> &ptr,
                                              const std::string &path, const bool verbose);
template void WriteVectorImage<VectorVolumeXF>(const itk::SmartPointer<VectorVolumeXF> &ptr,
                                               const std::string &path, const bool verbose);
template void WriteVectorImage<VectorVolumeI>(const itk::SmartPointer<VectorVolumeI> &ptr,
                                              const std::string &path, const bool verbose);
template void WriteVectorMagnitudeImage<VectorVolumeXF>(const VectorVolumeXF *ptr,
                                                        const std::string &   path,
                                                        const bool            verbose);
template void
              WriteVectorMagnitudeImage<VectorVolumeXF>(const itk::SmartPointer<VectorVolumeXF> &ptr,
                                          const std::string &path, const bool verbose);
template void WriteScaledVectorImage<VectorVolumeF>(const VectorVolumeF *img, const VolumeF *simg,
                                                    const std::string &path, const bool verbose);
template void WriteScaledVectorImage<VectorVolumeF>(const itk::SmartPointer<VectorVolumeF> &ptr,
                                                    const itk::SmartPointer<VolumeF> &      sptr,
                                                    const std::string &path, const bool verbose);
} // namespace QI
