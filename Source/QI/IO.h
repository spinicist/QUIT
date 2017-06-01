/*
 *  IO.h
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_IO_H

#include <string>
#include <sstream>

#include <Eigen/Core>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkDivideImageFilter.h"

#include "QI/Types.h"
#include "QI/Macro.h"

namespace QI {

template<typename TImg = QI::VolumeF>
auto ReadImage(const std::string &path) -> typename TImg::Pointer {
    typedef itk::ImageFileReader<TImg> TReader;
    typename TReader::Pointer file = TReader::New();
    file->SetFileName(path);
    file->Update();
    typename TImg::Pointer img = file->GetOutput();
    if (!img) {
        QI_EXCEPTION("Failed to read file: " << path);
    }
    img->DisconnectPipeline();
    return img;
}

template<typename TImg = QI::VolumeF>
auto ReadMagnitudeImage(const std::string &path) -> typename TImg::Pointer {
    typedef itk::Image<std::complex<typename TImg::PixelType>, TImg::ImageDimension> TComplex;
    auto x_img = QI::ReadImage<TComplex>(path);
    auto magFilter = itk::ComplexToModulusImageFilter<TComplex, TImg>::New();
    magFilter->SetInput(x_img);
    magFilter->Update();
    auto img = magFilter->GetOutput();
    img->DisconnectPipeline();
    return img;
}

template<typename TImg>
void WriteImage(const TImg *ptr, const std::string &path) {
    typedef itk::ImageFileWriter<TImg> TWriter;
    typename TWriter::Pointer file = TWriter::New();
    file->SetFileName(path);
    file->SetInput(ptr);
    file->Update();
}

template<typename TImg>
void WriteImage(const itk::SmartPointer<TImg> ptr, const std::string &path) {
    WriteImage<TImg>(ptr.GetPointer(), path);
}

template<typename TImg>
void WriteMagnitudeImage(const TImg *ptr, const std::string &path) {
    typedef typename TImg::PixelType::value_type TReal;
    typedef itk::Image<TReal, TImg::ImageDimension> TRealImage;
    auto mag = itk::ComplexToModulusImageFilter<TImg, TRealImage>::New();
    mag->SetInput(ptr);
    mag->Update();
    WriteImage<TRealImage>(mag->GetOutput(), path);
}

template<typename TImg>
void WriteMagnitudeImage(const itk::SmartPointer<TImg> ptr, const std::string &path) {
    WriteMagnitudeImage<TImg>(ptr.GetPointer(), path);
}

template<typename TImg>
void WriteScaledImage(const TImg *img, const QI::VolumeF *simg, const std::string &path) {
    auto scaleFilter = itk::DivideImageFilter<TImg, QI::VolumeF, TImg>::New();
    scaleFilter->SetInput1(img);
    scaleFilter->SetInput2(simg);
    scaleFilter->Update();
    WriteImage(scaleFilter->GetOutput(), path);
}

template<typename TImg>
void WriteScaledImage(const itk::SmartPointer<TImg> &ptr, const itk::SmartPointer<QI::VolumeF> &sptr, const std::string &path) {
    WriteScaledImage<TImg>(ptr.GetPointer(), sptr.GetPointer(), path);
}

template<typename TPixel>
auto ReadVectorImage(const std::string &path) -> typename itk::VectorImage<TPixel, 3>::Pointer {
    typedef itk::Image<TPixel, 4> TSeries;
    typedef itk::VectorImage<TPixel, 3> TVector;
    typedef itk::ImageToVectorFilter<TSeries> TToVector;
    
    auto img = ReadImage<TSeries>(path);
    auto convert = TToVector::New();
    convert->SetInput(img);
    convert->Update();
    typename TVector::Pointer vols = convert->GetOutput();
    vols->DisconnectPipeline();
    return vols;
}

template<typename TVImg>
void WriteVectorImage(const TVImg *img, const std::string &path) {
    typedef itk::VectorToImageFilter<TVImg> TToSeries;

    auto convert = TToSeries::New();
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

template<typename T> bool Read(const std::string &s, T &val) {
    std::istringstream stream(s);
    if (!(stream >> val)) {
        QI_EXCEPTION("Failed to parse input: " << s);
    }
    return true;
}

template<typename T> bool Read(std::istream &in, T &val) {
    std::string line;
    // Ignore comment lines. Use shell script convention
    while (in.peek() == '#') {
        if (!std::getline(in, line))
            QI_EXCEPTION("Failed to read input.");
    }
    if (!std::getline(in, line)) {
        QI_EXCEPTION("Failed to read input. Last line was: " << line);
    }
    return Read(line, val);
}

template<typename Scalar>
void ReadArray(const std::string &s, Eigen::Array<Scalar, Eigen::Dynamic, 1> &array) {
    std::istringstream stream(s);
    std::vector<Scalar> vals;

    Scalar temp;
    while (stream >> temp) {
        vals.push_back(temp);
    }

    array = Eigen::Array<Scalar, Eigen::Dynamic, 1>(vals.size());
    for (int i = 0; i < vals.size(); i++) {
        array[i] = vals[i];
    }
}

template<typename Scalar>
void ReadArray(std::istream &in, Eigen::Array<Scalar, Eigen::Dynamic, 1> &array) {
    std::string line;
    // Ignore comment lines. Use shell script convention
    while (in.peek() == '#') {
        if (!std::getline(in, line))
            QI_EXCEPTION("Failed to read input.");
    }
    if (!std::getline(in, line)) {
        QI_EXCEPTION("Failed to read input.");
    }
    ReadArray(line, array);
}

template<typename Scalar>
void ReadArray(const std::string &s, Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> &array) {
    std::istringstream stream(s);
    std::vector<Scalar> vals;

    Scalar temp;
    while (stream >> temp) {
        vals.push_back(temp);
    }

    array = Eigen::Array<Scalar, Eigen::Dynamic, 1>(vals.size());
    for (int i = 0; i < vals.size(); i++) {
        array[i] = vals[i];
    }
}

template<typename Scalar>
void ReadArray(std::istream &in, Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> &array) {
    std::vector<Eigen::Array<Scalar, Eigen::Dynamic, 1>> rows;
    std::string line;
    
    Eigen::Array<Scalar, Eigen::Dynamic, 1> row;
    while (std::getline(in, line) && line != "") {
        ReadArray(line, row);
        rows.push_back(row);
    }
    
    array = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>(rows.size(), rows.at(0).size());
    for (int i = 0; i < rows.size(); i++) {
        if (rows.at(i).size() == array.cols()) {
            array.row(i) = rows.at(i);
        } else {
            QI_EXCEPTION("Inconsistent row sizes in input");
        }
    }
}

} // End namespace QUIT

#endif // QUIT_IO_H