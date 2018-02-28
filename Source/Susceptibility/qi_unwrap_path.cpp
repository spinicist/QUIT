/*
 *  qi_unwrap_path.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  This is an implementation of the algorithm found in:
 *  Abdul-Rahman et al, Fast and robust three-dimensional best path phase unwrapping algorithm,
 *  http://ao.osa.org/abstract.cfm?URI=ao-46-26-6623
 *  Abdul-Rahman et al, Robust three-dimensional best-path phase-unwrapping algorithm that 
 *  avoids singularity loops, http://ao.osa.org/abstract.cfm?URI=ao-48-23-4582
 */

#include <iostream>
#include <memory>
#include <list>
#include <list>

#include "itkImageToImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkExtractImageFilter.h"
#include "itkTileImageFilter.h"

#include "ImageTypes.h"
#include "Util.h"
#include "ImageIO.h"
#include "Args.h"

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
    PhaseReliabilityFilter() {
        this->SetNumberOfRequiredInputs(1);
        this->SetNumberOfRequiredOutputs(1);
        this->SetNthOutput(0, this->MakeOutput(0));
    }
    ~PhaseReliabilityFilter() {}

    float wrap(float voxel_value) {
        if (voxel_value > M_PI) {
            return voxel_value - (2*M_PI);
        } else if (voxel_value < -M_PI) {
            return voxel_value + (2*M_PI);
        } else {
            return voxel_value;
        }
    }

    void ThreadedGenerateData(const TRegion &region, ThreadIdType threadId) ITK_OVERRIDE {
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
            for (auto j = 0; j < back.size(); j++) {
                const float d = wrap(inputIter.GetPixel(back[j]) - phase) -
                                wrap(phase - inputIter.GetPixel(fwrd[j]));
                reliability += d*d;
            }
            outputIter.Set(reliability);
            ++inputIter; ++outputIter;
        }
    }

private:
    PhaseReliabilityFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};


class UnwrapPathPhaseFilter : public ImageToImageFilter<QI::VolumeF, QI::VolumeF> {
public:
    typedef QI::VolumeF     TImage;
    typedef typename TImage::RegionType        RegionType;
    typedef UnwrapPathPhaseFilter              Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef SmartPointer<Self>                 Pointer;

    itkNewMacro(Self);
    itkTypeMacro(Self, Superclass);

    void SetReliability(const TImage *img) { this->SetNthInput(1, const_cast<TImage*>(img)); }
    void GenerateOutputInformation() ITK_OVERRIDE {
        Superclass::GenerateOutputInformation();
        auto op = this->GetOutput();
        op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
        op->Allocate();
    }

protected:
    UnwrapPathPhaseFilter() {
        this->SetNumberOfRequiredInputs(2);
        this->SetNumberOfRequiredOutputs(1);
        this->SetNthOutput(0, this->MakeOutput(0));
    }
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

    int find_wrap(float phase1, float phase2) {
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
    void merge_groups(std::shared_ptr<std::list<Voxel *>> group1,
                      std::shared_ptr<std::list<Voxel *>> group2,
                      int wrap_delta) {
        for (auto v : *group2) {
            v->group = group1;
            v->wraps += wrap_delta;
        }
        group1->splice(group1->end(), *group2);
    }

    void GenerateData() ITK_OVERRIDE {
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

private:
    UnwrapPathPhaseFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

} // End namespace itk

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Path-based phase unwrapping\n"
                                "See Abdul-Rahman et al. Fast and Robust three-dimensional best path phase unwrapping algorithm\n"
                                "http://ao.osa.org/abstract.cfm?URI=ao-46-26-6623\n"
                                "http://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "PHASE", "Wrapped phase image");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPUT PREFIX", "Change output prefix (default input filename)", {'o', "out"});
    args::ValueFlag<std::string> maskarg(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    QI::ParseArgs(parser, argc, argv);
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(threads.Get());

    if (verbose) std::cout << "Reading phase file: " << QI::CheckPos(input_path) << std::endl;
    auto inFile = QI::ReadImage<QI::SeriesF>(input_path.Get());
    if (verbose && maskarg) std::cout << "Reading mask file: " << maskarg.Get() << std::endl;
    auto maskFile = maskarg ? QI::ReadImage(maskarg.Get()) : ITK_NULLPTR;
    
    typedef itk::ExtractImageFilter<QI::SeriesF, QI::VolumeF> TExtract;
    typedef itk::TileImageFilter<QI::VolumeF, QI::SeriesF>    TTile;
    
    auto region = inFile->GetLargestPossibleRegion();
    const size_t nvols = region.GetSize()[3]; // Save for the loop
    region.GetModifiableSize()[3] = 0;
    itk::FixedArray<unsigned int, 4> layout;
    layout[0] = layout[1] = layout[2] = 1;
    layout[3] = nvols;

    auto extract = TExtract::New();
    auto tile   = TTile::New();
    extract->SetInput(inFile);
    extract->SetDirectionCollapseToSubmatrix();
    tile->SetLayout(layout);

    auto reliabilityFilter = itk::PhaseReliabilityFilter::New();
    auto unwrapFilter = itk::UnwrapPathPhaseFilter::New();
    for (int i = 0; i < nvols; i++) {
        region.GetModifiableIndex()[3] = i;
        if (verbose) std::cout << "Processing volume " << i << std::endl;
        extract->SetExtractionRegion(region);
        extract->Update();
        if (verbose) std::cout << "Calculating reliabilty" << std::endl;
        reliabilityFilter->SetInput(extract->GetOutput());
        reliabilityFilter->Update();
        if (verbose) std::cout << "Unwrapping phase" << std::endl;
        unwrapFilter->SetInput(extract->GetOutput());
        unwrapFilter->SetReliability(reliabilityFilter->GetOutput());
        unwrapFilter->Update();
        tile->SetInput(i, unwrapFilter->GetOutput());
        unwrapFilter->GetOutput()->DisconnectPipeline();
    }
    tile->Update();
    // Make sure output orientation info is correct
    auto dir = inFile->GetDirection();
    auto spc = inFile->GetSpacing();
    inFile = tile->GetOutput();
    inFile->SetDirection(dir);
    inFile->SetSpacing(spc);
    inFile->DisconnectPipeline();

    std::string outname = (outarg ? outarg.Get() : (QI::StripExt(input_path.Get())) + "_unwrapped" + QI::OutExt());
    if (verbose) std::cout << "Writing output: " << outname << std::endl;
    QI::WriteImage(inFile, outname);
    return EXIT_SUCCESS;
}
