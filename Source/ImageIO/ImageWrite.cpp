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

#include "ImageIO.h"
#include "Log.h"

namespace QI {

template <typename TImg>
void WriteImage(const TImg *ptr, const std::string &path, const bool verbose) {
    typedef itk::ImageFileWriter<TImg> TWriter;
    typename TWriter::Pointer          file = TWriter::New();
    file->SetFileName(path);
    file->SetInput(ptr);
    QI::Log(verbose, "Writing image: {}", path);
    file->Update();
}

template <typename TImg>
void WriteImage(const itk::SmartPointer<TImg> &ptr, const std::string &path, const bool verbose) {
    WriteImage<TImg>(ptr.GetPointer(), path, verbose);
}

template <typename TImg>
void WriteMagnitudeImage(const TImg *ptr, const std::string &path, const bool verbose) {
    typedef typename TImg::PixelType::value_type    TReal;
    typedef itk::Image<TReal, TImg::ImageDimension> TRealImage;
    auto mag = itk::ComplexToModulusImageFilter<TImg, TRealImage>::New();
    mag->SetInput(ptr);
    mag->Update();
    WriteImage<TRealImage>(mag->GetOutput(), path, verbose);
}

template <typename TImg>
void WriteMagnitudeImage(const itk::SmartPointer<TImg> &ptr, const std::string &path,
                         const bool verbose) {
    WriteMagnitudeImage<TImg>(ptr.GetPointer(), path, verbose);
}

template <typename TImg>
void WriteScaledImage(const TImg *img, const QI::VolumeF *simg, const std::string &path,
                      const bool verbose) {
    auto scaleFilter = itk::DivideImageFilter<TImg, QI::VolumeF, TImg>::New();
    scaleFilter->SetInput1(img);
    scaleFilter->SetInput2(simg);
    scaleFilter->Update();
    WriteImage(scaleFilter->GetOutput(), path, verbose);
}

template <typename TImg>
void WriteScaledImage(const itk::SmartPointer<TImg> &       ptr,
                      const itk::SmartPointer<QI::VolumeF> &sptr, const std::string &path,
                      const bool verbose) {
    WriteScaledImage<TImg>(ptr.GetPointer(), sptr.GetPointer(), path, verbose);
}

template void WriteImage<VolumeF>(const VolumeF *ptr, const std::string &path, const bool verbose);
template void WriteImage<VolumeXF>(const VolumeXF *ptr, const std::string &path,
                                   const bool verbose);
template void WriteImage<VolumeD>(const VolumeD *ptr, const std::string &path, const bool verbose);
template void WriteImage<VolumeI>(const VolumeI *ptr, const std::string &path, const bool verbose);
template void WriteImage<VolumeUC>(const VolumeUC *ptr, const std::string &path,
                                   const bool verbose);
template void WriteImage<SeriesF>(const SeriesF *ptr, const std::string &path, const bool verbose);
template void WriteImage<SeriesD>(const SeriesD *ptr, const std::string &path, const bool verbose);
template void WriteImage<SeriesI>(const SeriesI *ptr, const std::string &path, const bool verbose);
template void WriteImage<SeriesXF>(const SeriesXF *ptr, const std::string &path,
                                   const bool verbose);
template void WriteImage<SeriesXD>(const SeriesXD *ptr, const std::string &path,
                                   const bool verbose);
template void WriteImage<VolumeF>(const itk::SmartPointer<VolumeF> &ptr, const std::string &path,
                                  const bool verbose);
template void WriteImage<VolumeD>(const itk::SmartPointer<VolumeD> &ptr, const std::string &path,
                                  const bool verbose);
template void WriteImage<VolumeI>(const itk::SmartPointer<VolumeI> &ptr, const std::string &path,
                                  const bool verbose);
template void WriteImage<SeriesF>(const itk::SmartPointer<SeriesF> &ptr, const std::string &path,
                                  const bool verbose);
template void WriteImage<SeriesD>(const itk::SmartPointer<SeriesD> &ptr, const std::string &path,
                                  const bool verbose);
template void WriteImage<SeriesXF>(const itk::SmartPointer<SeriesXF> &ptr, const std::string &path,
                                   const bool verbose);
template void WriteImage<SeriesXD>(const itk::SmartPointer<SeriesXD> &ptr, const std::string &path,
                                   const bool verbose);
template void WriteScaledImage<VolumeF>(const VolumeF *img, const VolumeF *simg,
                                        const std::string &path, const bool verbose);
template void WriteScaledImage<VolumeF>(const itk::SmartPointer<VolumeF> &ptr,
                                        const itk::SmartPointer<VolumeF> &sptr,
                                        const std::string &path, const bool verbose);
template void WriteMagnitudeImage<VolumeXF>(const VolumeXF *ptr, const std::string &path,
                                            const bool verbose);
template void WriteMagnitudeImage<VolumeXF>(const itk::SmartPointer<VolumeXF> &ptr,
                                            const std::string &path, const bool verbose);
template void WriteMagnitudeImage<SeriesXF>(const SeriesXF *ptr, const std::string &path,
                                            const bool verbose);
template void WriteMagnitudeImage<SeriesXF>(const itk::SmartPointer<SeriesXF> &ptr,
                                            const std::string &path, const bool verbose);

} // namespace QI
