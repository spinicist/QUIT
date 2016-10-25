/*
 *  MeanImageFilter.hxx
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef MEANIMAGEFILTER_HXX
#define MEANIMAGEFILTER_HXX

namespace itk {

template<typename TInput>
MeanImageFilter<TInput>::MeanImageFilter() {}

template<typename TInput>
void MeanImageFilter<TInput>::SetNumberOfImages(const size_t n) {
    //std::cout <<  __PRETTY_FUNCTION__ << endl;
    this->SetNumberOfRequiredInputs(n);
    this->SetNumberOfRequiredOutputs(1);
    this->SetNthOutput(0, this->MakeOutput(0));
}

template<typename TInput>
void MeanImageFilter<TInput>::GenerateOutputInformation() {
    //std::cout << __PRETTY_FUNCTION__ << std::endl;
    typename Superclass::OutputImagePointer outputPtr = this->GetOutput();
    typename Superclass::InputImageConstPointer inputPtr  = this->GetInput();
    if ( !outputPtr || !inputPtr ) {
        return;
    }

    outputPtr->SetNumberOfComponentsPerPixel(inputPtr->GetNumberOfComponentsPerPixel());
    outputPtr->SetRegions(inputPtr->GetLargestPossibleRegion());
    outputPtr->SetSpacing(inputPtr->GetSpacing());
    outputPtr->SetOrigin(inputPtr->GetOrigin());
    outputPtr->SetDirection(inputPtr->GetDirection());
    outputPtr->Allocate();
    //std::cout << "END " << __PRETTY_FUNCTION__ << std::endl;
}

template<typename TInput>
void MeanImageFilter<TInput>::ThreadedGenerateData(const TRegion &region, ThreadIdType threadId) {
    //std::cout << __PRETTY_FUNCTION__ << std::endl;
    size_t N = this->GetNumberOfRequiredInputs();
    std::vector<ImageRegionConstIterator<TInput>> iters(N);
    for (size_t i = 0; i < N; i++) {
        iters[i] = ImageRegionConstIterator<TInput>(this->GetInput(i), region);
    }
    ImageRegionIterator<TOutput> outputIter(this->GetOutput(), region);
    //TimeProbe clock;
    //clock.Start();
    while(!outputIter.IsAtEnd()) {
        double accumulator = 0.;
        for (auto& it : iters) {
            accumulator += it.Get() / N;
            ++it;
        }
        outputIter.Set(accumulator);
        ++outputIter;
    }
    //clock.Stop();
    //std::cout << "END " << __PRETTY_FUNCTION__ << std::endl;
}

} // End namespace itk

#endif // Define MEANIMAGEFILTER_HXX
