/*
 *  qi_zshim.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Core>

#include "Util.h"
#include "ImageIO.h"
#include "Args.h"
#include "itkUnaryFunctorImageFilter.h"

struct SoSFunctor {
    int m_outsize = 0, m_zshims = 0;

    SoSFunctor() = default;
    SoSFunctor(const int isz, const int z) :
        m_outsize(isz / z),
        m_zshims(z)
    {};
    SoSFunctor &operator=(const SoSFunctor &other) {
        m_outsize = other.m_outsize;
        m_zshims = other.m_zshims;
        return *this;
    }

    ~SoSFunctor() {};
    bool operator!=(const SoSFunctor &) const { return true; }
    bool operator==(const SoSFunctor &other) const { return !(*this != other); }
    inline itk::VariableLengthVector<float> operator()(const itk::VariableLengthVector<float> &vec) const {
        Eigen::Map<const Eigen::MatrixXf> input(vec.GetDataPointer(), m_zshims, m_outsize);
        const Eigen::VectorXf output = input.colwise().norm();
        itk::VariableLengthVector<float> out_vec(m_outsize);
        for (int i = 0; i < m_outsize; i++) {
            out_vec[i] = output[i];
        }
        return out_vec;
    }
};

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Combines Z-Shimmed data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "ZSHIM_FILE", "Input Z-Shimmed file");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Output progress messages", {'v', "verbose"});
    args::Flag     debug(parser, "DEBUG", "Output debug messages", {'d', "debug"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    args::ValueFlag<int> zshims(parser, "ZSHIMS", "Number of Z-Shims (default 8)", {'z', "zshims"}, 8);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    auto input = QI::ReadVectorImage(QI::CheckPos(input_path), verbose);
    const std::string outPrefix = outarg ? outarg.Get() : QI::Basename(input_path.Get());
    SoSFunctor sos(input->GetNumberOfComponentsPerPixel(), zshims.Get());
    auto sos_filter = itk::UnaryFunctorImageFilter<QI::VectorVolumeF, QI::VectorVolumeF, SoSFunctor>::New();
    sos_filter->SetFunctor(sos);
    sos_filter->SetInput(input);
    QI_LOG(verbose, "Processing");
    sos_filter->Update();
    const std::string fname = outPrefix + "_zshim" + QI::OutExt();
    QI_LOG(verbose, "Writing file: " << fname);
    QI::WriteVectorImage(sos_filter->GetOutput(), fname);
    return EXIT_SUCCESS;
}
