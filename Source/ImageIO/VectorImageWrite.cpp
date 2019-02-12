/*
 *  ImageWrite.cpp
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
void WriteImage(const TVImg *img, const std::string &path, const bool verbose) {
    using TToSeries = itk::VectorToImageFilter<TVImg>;
    using TWriter   = itk::ImageFileWriter<typename TToSeries::TOutput>;

    typename TToSeries::Pointer convert = TToSeries::New();
    convert->SetInput(img);
    convert->Update();

    typename TWriter::Pointer file = TWriter::New();
    file->SetFileName(path);
    file->SetInput(convert->GetOutput());
    QI::Log(verbose, "Writing image: {}", path);
    file->Update();
}

template <typename TVImg>
void WriteImage(const itk::SmartPointer<TVImg> &ptr, const std::string &path, const bool verbose) {
    WriteImage(ptr.GetPointer(), path, verbose);
}

template <typename TVImg>
void WriteMagnitudeImage(const TVImg *img, const std::string &path, const bool verbose) {
    using TToSeries   = itk::VectorToImageFilter<TVImg>;
    using TPixel      = typename TVImg::InternalPixelType;
    using TReal       = typename TPixel::value_type;
    using TSeries     = itk::Image<TPixel, 4>;
    using TRealSeries = itk::Image<TReal, 4>;

    auto convert = TToSeries::New();
    convert->SetInput(img);

    auto mag = itk::ComplexToModulusImageFilter<TSeries, TRealSeries>::New();
    mag->SetInput(convert->GetOutput());
    mag->Update();

    using TWriter = itk::ImageFileWriter<TRealSeries>;
    auto file     = TWriter::New();
    file->SetFileName(path);
    file->SetInput(mag->GetOutput());
    QI::Log(verbose, "Writing magnitude image: {}", path);
    file->Update();
}

template <typename TVImg>
void WriteMagnitudeImage(const itk::SmartPointer<TVImg> &ptr, const std::string &path,
                         const bool verbose) {
    WriteMagnitudeImage(ptr.GetPointer(), path, verbose);
}

template <typename TVImg>
void WriteScaledImage(const TVImg *img, const QI::VolumeF *simg, const std::string &path,
                      const bool verbose) {
    auto scaleFilter = itk::DivideImageFilter<TVImg, QI::VolumeF, TVImg>::New();
    scaleFilter->SetInput1(img);
    scaleFilter->SetInput2(simg);
    scaleFilter->Update();
    WriteImage(scaleFilter->GetOutput(), path, verbose);
}

template <typename TVImg>
void WriteScaledImage(const itk::SmartPointer<TVImg> &      ptr,
                      const itk::SmartPointer<QI::VolumeF> &sptr, const std::string &path,
                      const bool verbose) {
    WriteScaledImage(ptr.GetPointer(), sptr.GetPointer(), path, verbose);
}

template void WriteImage<VectorVolumeF>(const VectorVolumeF *img, const std::string &path,
                                        const bool verbose);
template void WriteImage<VectorVolumeXF>(const VectorVolumeXF *img, const std::string &path,
                                         const bool verbose);
template void WriteImage<VectorVolumeF>(const itk::SmartPointer<VectorVolumeF> &ptr,
                                        const std::string &path, const bool verbose);
template void WriteImage<VectorVolumeXF>(const itk::SmartPointer<VectorVolumeXF> &ptr,
                                         const std::string &path, const bool verbose);
template void WriteImage<VectorVolumeI>(const itk::SmartPointer<VectorVolumeI> &ptr,
                                        const std::string &path, const bool verbose);
template void WriteMagnitudeImage<VectorVolumeXF>(const VectorVolumeXF *ptr,
                                                  const std::string &path, const bool verbose);
template void WriteMagnitudeImage<VectorVolumeXF>(const itk::SmartPointer<VectorVolumeXF> &ptr,
                                                  const std::string &path, const bool verbose);
template void WriteScaledImage<VectorVolumeF>(const VectorVolumeF *img, const VolumeF *simg,
                                              const std::string &path, const bool verbose);
template void WriteScaledImage<VectorVolumeF>(const itk::SmartPointer<VectorVolumeF> &ptr,
                                              const itk::SmartPointer<VolumeF> &      sptr,
                                              const std::string &path, const bool verbose);
} // namespace QI
