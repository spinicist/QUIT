/*
 *  qi_lineshape.cpp - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2018 Tobias Wood, Samuel Hurley, Erika Raven
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Args.h"
#include "JSON.h"
#include "Lineshape.h"
#include "Util.h"

int lineshape_main(int argc, char **argv) {
    Eigen::initParallel();

    args::ArgumentParser parser(
        "Calculates different lineshapes then saves them as a splineshape for fast qMT lookup.\n"
        "http://github.com/spinicist/QUIT");

    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Positional<std::string> out_path(parser, "OUTPUT", "Output lineshape JSON path");
    args::ValueFlag<std::string>  shape_arg(parser,
                                           "LINESHAPE",
                                           "Choose lineshape (Gauss/Lorentzian/SuperLorentzian)",
                                           {'l', "lineshape"},
                                           "G");
    args::ValueFlag<double>       T2b(
        parser, "T2 BOUND", "Choose nominal bound-pool T2 (default 10µs)", {'t', "T2b"}, 10e-6);
    args::ValueFlag<double> frq_count(
        parser, "FREQUENCY COUNT", "Number of frequencies (default 10)", {'n', "frq_count"}, 10);
    args::ValueFlag<double> frq_start(parser,
                                      "FREQUENCY START",
                                      "First saturation frequency (default 1000 Hz)",
                                      {'s', "frq_start"},
                                      1e3);
    args::ValueFlag<double> frq_spacing(parser,
                                        "FREQUENCY SPACING",
                                        "Spacing of frequencies (default 1000 Hz)",
                                        {'p', "frq_space"},
                                        1e3);

    QI::ParseArgs(parser, argc, argv, verbose);
    QI::Log(verbose, "Bound-pool T2: {}", T2b.Get());
    auto frqs =
        Eigen::ArrayXd::LinSpaced(frq_count.Get(),
                                  frq_start.Get(),
                                  frq_start.Get() + frq_spacing.Get() * (frq_count.Get() - 1));
    QI::Log(verbose, "Frequency count: {}", frqs.rows());
    Eigen::ArrayXd values;
    if (shape_arg.Get() == "Gaussian") {
        values = QI::Gaussian(frqs, T2b.Get());
    } else if (shape_arg.Get() == "Lorentzian") {
        values = QI::Lorentzian(frqs, T2b.Get());
    } else if (shape_arg.Get() == "SuperLorentzian") {
        values = QI::SuperLorentzian(frqs, T2b.Get());
    } else {
        QI::Fail("Unknown lineshape: {}", shape_arg.Get());
    }
    const auto lineshape =
        QI::InterpLineshape(frq_start.Get(), frq_spacing.Get(), frq_count.Get(), values, T2b.Get());
    json doc{{"lineshape", lineshape}};
    QI::WriteJSON(QI::CheckPos(out_path), doc);
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
