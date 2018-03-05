/*
 *  PathUnwrapFilter.cpp
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

#include "PathUnwrapFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"

namespace itk {

void UnwrapPathPhaseFilter::SetReliability(const TImage *img) { this->SetNthInput(1, const_cast<TImage*>(img)); }
void UnwrapPathPhaseFilter::GenerateOutputInformation() {
    Superclass::GenerateOutputInformation();
    auto op = this->GetOutput();
    op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    op->Allocate();
}

UnwrapPathPhaseFilter::UnwrapPathPhaseFilter() {
    this->SetNumberOfRequiredInputs(2);
    this->SetNumberOfRequiredOutputs(1);
    this->SetNthOutput(0, this->MakeOutput(0));
}

int UnwrapPathPhaseFilter::find_wrap(float phase1, float phase2) {
    const float difference = phase1 - phase2;
    if (difference > M_PI) {
        return -1;
    } else if (difference < -M_PI) {
        return 1;
    } else {
        return 0;
    }
}

// Do not point shared_ptrs as references, otherwise group2 gets deallocated during loop
void UnwrapPathPhaseFilter::merge_groups(std::shared_ptr<std::list<Voxel *>> group1,
                    std::shared_ptr<std::list<Voxel *>> group2,
                    int wrap_delta) {
    for (auto v : *group2) {
        v->group = group1;
        v->wraps += wrap_delta;
    }
    group1->splice(group1->end(), *group2);
}

void UnwrapPathPhaseFilter::GenerateData() {
    //std::cout <<  __PRETTY_FUNCTION__ << std::endl;
    const auto region = this->GetInput()->GetLargestPossibleRegion();
    const auto volume_size = region.GetSize()[0] * region.GetSize()[1] * region.GetSize()[2];
    int volume_width = region.GetSize()[0];
    int volume_height = region.GetSize()[1];
    int volume_depth = region.GetSize()[2];

    std::vector<Voxel> voxels;
    voxels.reserve(volume_size);
    ImageRegionConstIterator<TImage> phaseIter(this->GetInput(0), region);
    ImageRegionConstIterator<TImage> relIter(this->GetInput(1), region);
    while (!phaseIter.IsAtEnd()) {
        voxels.emplace_back(Voxel{phaseIter.Get(), relIter.Get()});
        voxels.back().group->push_front(&voxels.back());
        ++phaseIter; ++relIter;
    }
    std::vector<Edge> edges;
    edges.reserve(3 * volume_size);
    std::vector<Voxel>::iterator v_it = voxels.begin();
    for (int n=0; n < volume_depth; n++) {
        for (int i = 0; i < volume_height; i++) {
            for (int j = 0; j < volume_width - 1; j++) {
                edges.emplace_back(Edge{v_it->reliability + (v_it + 1)->reliability,
                                    &*v_it, &*(v_it + 1),
                                    find_wrap(v_it->phase, (v_it + 1)->phase)});
                v_it++;
            }
        v_it++;
        }
    }

    v_it = voxels.begin();
    for (int n=0; n < volume_depth; n++) {
        for (int i = 0; i < volume_height - 1; i++) {
            for (int j = 0; j < volume_width; j++) {
                edges.emplace_back(Edge{v_it->reliability + (v_it + volume_width)->reliability,
                                    &*v_it, &*(v_it + volume_width),
                                    find_wrap(v_it->phase, (v_it + volume_width)->phase)});
                v_it++;
            }
        }
        v_it+=volume_width;
    }

    v_it = voxels.begin();
    for (int n=0; n < volume_depth - 1; n++) {
        for (int i = 0; i < volume_height; i++) {
            for (int j = 0; j < volume_width; j++) {
                edges.emplace_back(Edge{v_it->reliability + (v_it + volume_width*volume_height)->reliability,
                                    &*v_it, &*(v_it + (volume_width*volume_height)),
                                    find_wrap(v_it->phase, (v_it + volume_width*volume_height)->phase)});
                v_it++;
            }
        }
    }

    std::stable_sort(edges.begin(), edges.end(), [](Edge a, Edge b){ return a.reliability < b.reliability; });
    for (auto &edge : edges) {
        Voxel *voxel1 = edge.voxel1;
        Voxel *voxel2 = edge.voxel2;
        if (voxel2->group != voxel1->group) {
            if (voxel1->group->size() > voxel2->group->size()) {
                merge_groups(voxel1->group, voxel2->group, voxel1->wraps - edge.wrap - voxel2->wraps);
            } else {
                merge_groups(voxel2->group, voxel1->group, voxel2->wraps + edge.wrap - voxel1->wraps);
            }
        }
    }

    // Unwrap voxels and reassemble into image
    ImageRegionIterator<TImage> outputIter(this->GetOutput(), region);
    v_it = voxels.begin();
    while (!outputIter.IsAtEnd()) {
        outputIter.Set(v_it->phase + 2*M_PI*v_it->wraps);
        ++v_it; ++outputIter;
    }
}

} // End namespace itk