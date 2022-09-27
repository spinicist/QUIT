/*
 *  PathUnwrapFilter.h
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

#ifndef PATH_UNWRAP_FILTER_H
#define PATH_UNWRAP_FILTER_H

#include <memory>
#include "itkImageToImageFilter.h"
#include "ImageTypes.h"

namespace itk {

class UnwrapPathPhaseFilter : public ImageToImageFilter<QI::VolumeF, QI::VolumeF> {
public:
    typedef QI::VolumeF     TImage;
    typedef typename TImage::RegionType        RegionType;
    typedef UnwrapPathPhaseFilter              Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef SmartPointer<Self>                 Pointer;

    itkNewMacro(Self)
    itkTypeMacro(Self, Superclass);

    void SetReliability(const TImage *img);
    void GenerateOutputInformation() ITK_OVERRIDE;

protected:
    UnwrapPathPhaseFilter();
    ~UnwrapPathPhaseFilter() {}

    struct Voxel {
        float phase;
        float reliability;
        int wraps;
        std::shared_ptr<std::list<Voxel *>> group;

        Voxel() {}
        Voxel(float ph, float rel) {
            phase = ph; reliability = rel;
            wraps = 0;
            group = std::make_shared<std::list<Voxel *>>();
            // group->push_front(this);
            // std::cout << "New Voxel " << this << " / " << group->front() << std::endl;
        }
    };

    struct Edge {
        float reliability;			//reliabilty of the edge and its equal to the sum of the reliability of the two voxels that it conntects
        Voxel *voxel1;		//pointer to the first voxel
        Voxel *voxel2;		//pointer to the second voxel
        int wrap;			//No. of 2*pi to add to the voxel to unwrap it 
    }; 

    int find_wrap(float phase1, float phase2);

    // Do not point shared_ptrs as references, otherwise group2 gets deallocated during loop
    void merge_groups(std::shared_ptr<std::list<Voxel *>> group1,
                      std::shared_ptr<std::list<Voxel *>> group2,
                      int wrap_delta);

    void GenerateData() ITK_OVERRIDE;

private:
    UnwrapPathPhaseFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

} // End namespace itk

#endif // PATH_UNWRAP_FILTER_H