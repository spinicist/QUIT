/*
 *  ReliabilityFilter.cpp
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

#include "ReliabilityFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"

namespace itk {

PhaseReliabilityFilter::PhaseReliabilityFilter() {
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(1);
    this->SetNthOutput(0, this->MakeOutput(0));
}

float PhaseReliabilityFilter::wrap(float voxel_value) {
    if (voxel_value > M_PI) {
        return voxel_value - (2*M_PI);
    } else if (voxel_value < -M_PI) {
        return voxel_value + (2*M_PI);
    } else {
        return voxel_value;
    }
}

void PhaseReliabilityFilter::DynamicThreadedGenerateData(const TRegion &region) {
    using TIter = ImageRegionIterator<TImage>;
    using TCNIter = ConstNeighborhoodIterator<TImage>;
    using TOffset = TCNIter::OffsetType;
    TCNIter::RadiusType radius; radius.Fill(1);
    TCNIter inputIter(radius, this->GetInput(0), region);
    TIter   outputIter(this->GetOutput(), region);
    std::vector<TOffset> back = { {{-1, 0, 0}}, {{ 0,-1, 0}}, {{ 0, 0,-1}},
                                    {{-1,-1, 0}}, {{ 1,-1, 0}}, {{-1,-1,-1}},
                                    {{ 0,-1,-1}}, {{ 1,-1,-1}}, {{-1, 0,-1}},
                                    {{-1, 1,-1}}, {{ 1, 0,-1}}, {{ 0, 1,-1}},
                                    {{ 1, 1,-1}} };
    std::vector<TOffset> fwrd = { {{ 1, 0, 0}}, {{ 0, 1, 0}}, {{ 0, 0, 1}},
                                    {{ 1, 1, 0}}, {{-1, 1, 0}}, {{ 1, 1, 1}},
                                    {{ 0, 1, 1}}, {{-1, 1, 1}}, {{ 1, 0, 1}},
                                    {{ 1,-1, 1}}, {{-1, 0, 1}}, {{ 0,-1, 1}},
                                    {{-1,-1, 1}} };

    inputIter.GoToBegin();
    outputIter.GoToBegin();

    while (!inputIter.IsAtEnd()) {
        const float phase = inputIter.GetCenterPixel();
        float reliability = 0;
        for (size_t j = 0; j < back.size(); j++) {
            const float d = wrap(inputIter.GetPixel(back[j]) - phase) -
                            wrap(phase - inputIter.GetPixel(fwrd[j]));
            reliability += d*d;
        }
        outputIter.Set(reliability);
        ++inputIter; ++outputIter;
    }
}

} // End namespace itk