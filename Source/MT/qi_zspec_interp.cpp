/*
 *  qi_cestasym.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Core>

#include "Util.h"
#include "Args.h"
#include "ImageIO.h"
#include "JSON.h"
#include "Spline.h"
#include "itkBinaryFunctorImageFilter.h"

    
struct InterpFunctor {
    Eigen::ArrayXd m_ifrqs, m_ofrqs;

    InterpFunctor() = default;
    InterpFunctor(const Eigen::ArrayXd &ifrq, const Eigen::ArrayXd &ofrq) :
        m_ifrqs(ifrq),
        m_ofrqs(ofrq)
    {};

    ~InterpFunctor() {};
    bool operator!=(const InterpFunctor &) const { return true; }
    bool operator==(const InterpFunctor &other) const { return !(*this != other); }
    inline itk::VariableLengthVector<float> operator()(const itk::VariableLengthVector<float> &vec, const float &f0) const {
        const Eigen::Map<const Eigen::ArrayXf> zdata(vec.GetDataPointer(), vec.Size());
        QI::SplineInterpolator zspec(m_ifrqs, zdata.cast<double>());
        itk::VariableLengthVector<float> output(m_ofrqs.rows());
        for (int f = 0; f < m_ofrqs.rows(); f++) {
            const double frq = m_ofrqs[f] + f0;
            output[f] = zspec(frq);
        }
        return output;
    }
};

int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Interpolates Z-spectrums using Cubic Splines\nhttp://github.com/spinicist/QUIT");
    
    args::Positional<std::string> input_path(parser, "INPUT", "Input Z-spectrum file");
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPUT", "Change ouput filename (default is input_interp)", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> f0(parser, "OFF RESONANCE", "Specify off-resonance frequency (units must match input)", {'f', "f0"});
    args::ValueFlag<std::string> subregion(parser, "REGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI_LOG(verbose, "Opening file: " << QI::CheckPos(input_path));
    auto data = QI::ReadVectorImage<float>(QI::CheckPos(input_path), verbose);

    QI_LOG(verbose, "Reading input frequences");
    rapidjson::Document json = QI::ReadJSON(std::cin);
    auto i_frqs = QI::ArrayFromJSON(json["input_freqs"]);
    auto o_frqs = QI::ArrayFromJSON(json["output_freqs"]);
    InterpFunctor functor(i_frqs, o_frqs);
    auto filter = itk::BinaryFunctorImageFilter<QI::VectorVolumeF,
                                                QI::VolumeF,
                                                QI::VectorVolumeF,
                                                InterpFunctor>::New();
    filter->SetFunctor(functor);
    filter->SetInput1(data);
    filter->SetInput2(QI::ReadImage(f0.Get(), verbose));
    QI_LOG(verbose, "Processing");
    filter->Update();
    std::string outname = outarg ? outarg.Get() : QI::StripExt(input_path.Get()) + "_interp" + QI::OutExt();
    QI_LOG(verbose, "Writing output: " << outname);
    QI::WriteVectorImage(filter->GetOutput(), outname);
    QI_LOG(verbose, "Finished.");
    return EXIT_SUCCESS;
}
