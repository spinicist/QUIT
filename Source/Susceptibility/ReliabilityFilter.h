/*
 *  ReliabilityFilter.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  This is part of an implementation of the algorithm found in:
 *  Abdul-Rahman et al, Fast and robust three-dimensional best path phase unwrapping algorithm,
 *  http://ao.osa.org/abstract.cfm?URI=ao-46-26-6623
 *  Abdul-Rahman et al, Robust three-dimensional best-path phase-unwrapping algorithm that 
 *  avoids singularity loops, http://ao.osa.org/abstract.cfm?URI=ao-48-23-4582
 */

#ifndef RELIABILITY_FILTER_H
#define RELIABILITY_FILTER_H

#include "itkImageToImageFilter.h"
#include "ImageTypes.h"

namespace itk {

class PhaseReliabilityFilter : public ImageToImageFilter<QI::VolumeF, QI::VolumeF> {
public:
    using TImage     = QI::VolumeF;
    using TRegion    = TImage::RegionType;
    using Self       = PhaseReliabilityFilter;
    using Superclass = ImageToImageFilter<TImage, TImage>;
    using Pointer    = SmartPointer<Self>;

    itkNewMacro(Self);
    itkTypeMacro(Self, Superclass);

protected:
    PhaseReliabilityFilter();
    ~PhaseReliabilityFilter() {}

    float wrap(float voxel_value);
    void DynamicThreadedGenerateData(const TRegion &region) ITK_OVERRIDE;

private:
    PhaseReliabilityFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

} // End namespace ITK

#endif // RELIABILITY_FILTER_H