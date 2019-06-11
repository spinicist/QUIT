/*
 * ImageTypes.h
 *
 * Copyright (c) 2015, 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_IMAGE_TYPES
#define QUIT_IMAGE_TYPES

#include "itkImage.h"
#include "itkVectorImage.h"

namespace QI {

typedef itk::Image<unsigned char, 3> VolumeUC;
typedef itk::Image<unsigned int, 3>  VolumeUI;
typedef itk::Image<int, 3>           VolumeI;
typedef itk::Image<int, 4>           SeriesI;
typedef itk::VectorImage<int, 3>     VectorVolumeI;

typedef itk::Image<float, 3>                     VolumeF;
typedef itk::Image<float, 4>                     SeriesF;
typedef itk::VectorImage<float, 3>               VectorVolumeF;
typedef itk::Image<std::complex<float>, 3>       VolumeXF;
typedef itk::Image<std::complex<float>, 4>       SeriesXF;
typedef itk::VectorImage<std::complex<float>, 3> VectorVolumeXF;

typedef itk::Image<double, 3>                     VolumeD;
typedef itk::Image<double, 4>                     SeriesD;
typedef itk::VectorImage<double, 3>               VectorVolumeD;
typedef itk::Image<std::complex<double>, 3>       VolumeXD;
typedef itk::Image<std::complex<double>, 4>       SeriesXD;
typedef itk::VectorImage<std::complex<double>, 3> VectorVolumeXD;

template <typename TNew = VolumeF, typename TRef>
auto NewImageLike(const itk::SmartPointer<TRef> &ref, int const &ncomp = 1) ->
    typename TNew::Pointer {
    auto nimg = TNew::New();
    nimg->CopyInformation(ref);
    nimg->SetRegions(ref->GetBufferedRegion());
    if (ncomp > 1) {
        nimg->SetNumberOfComponentsPerPixel(ncomp);
    }
    nimg->Allocate(true);
    return nimg;
}

} // End namespace QI

#endif // define QUIT_IMAGE_TYPES
