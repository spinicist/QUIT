/*
 *  qi_mtsat.cpp
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>

#include "Args.h"
#include "ImageIO.h"
#include "MTSatModel.h"
#include "ModelFunc.h"
#include "SimulateModel.h"
#include "Util.h"

/*
 * Main
 */
int mtsat_main(args::Subparser &parser) {
    args::Positional<std::string> pdw_path(parser, "PDw", "Input PD-weighted file");
    args::Positional<std::string> t1w_path(parser, "T1w", "Input T1-weighted file");
    args::Positional<std::string> mtw_path(parser, "MTw", "Input MT-weighted file");
    QI_COMMON_ARGS;
    args::ValueFlag<std::string>  b1_path(parser, "B1", "Path to B1 map", {'b', "B1"});
    args::ValueFlag<double>       C(
        parser, "C", "Correction factor for delta (default 0.4)", {'C', "C"}, 0.4);
    args::Flag                    smallangle(
        parser, "smallangle", "Use small flip angle approx for R1 and PD calculation", {'s', "smallangle"});
    args::ValueFlag<double>       delta_max(
        parser, "delta_max", "Values of delta (MTsat) above this are clamped, in % (default 10%)", {'d', "delta-max"}, 10);
    args::ValueFlag<double>       r1_max(
        parser, "r1_max", "Values of R1 above this are clamped, in 1/s (default 10 1/s)", {'r', "r1-max"}, 10);
    parser.Parse();
    
    QI::Log(verbose, "Reading sequence parameters");
    json              input = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    QI::MTSatSequence seq(input["MTSat"]);

    MTSatModel model{{}, seq, C.Get(), smallangle, delta_max.Get(), r1_max.Get()};
    if (simulate) {
        QI::SimulateModel<MTSatModel, true>(
            input,
            model,
            {b1_path.Get()},
            {QI::CheckPos(pdw_path), QI::CheckPos(t1w_path), QI::CheckPos(mtw_path)},
            mask.Get(),
            verbose,
            simulate.Get(),
            threads.Get(),
            subregion.Get());
    } else {
        auto pdw_img = QI::ReadImage(QI::CheckPos(pdw_path), verbose);
        auto t1w_img = QI::ReadImage(QI::CheckPos(t1w_path), verbose);
        auto mtw_img = QI::ReadImage(QI::CheckPos(mtw_path), verbose);

        QI::VolumeF::Pointer B1_img = nullptr, mask_img = nullptr;
        if (b1_path) {
            B1_img = QI::ReadImage(b1_path.Get(), verbose);
        }
        if (mask) {
            mask_img = QI::ReadImage(mask.Get(), verbose);
        }

        QI::Info(verbose, "Allocating output memory");
        auto A_img  = QI::NewImageLike<QI::VolumeF>(pdw_img);
        auto R1_img = QI::NewImageLike<QI::VolumeF>(pdw_img);
        auto d_img  = QI::NewImageLike<QI::VolumeF>(pdw_img);

        std::array<QI::VolumeF::Pointer, 3> outs{A_img, R1_img, d_img};

        QI::ModelFunc(model, {pdw_img, t1w_img, mtw_img}, {B1_img}, mask_img, threads.Get(), outs);
        QI::Info(verbose, "Finished");
        QI::WriteImage(A_img, prefix.Get() + "MTSat_PD" + QI::OutExt(), verbose);
        QI::WriteImage(R1_img, prefix.Get() + "MTSat_R1" + QI::OutExt(), verbose);
        QI::WriteImage(d_img, prefix.Get() + "MTSat_delta" + QI::OutExt(), verbose);
    }
    return EXIT_SUCCESS;
}
