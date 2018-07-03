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

#include <iostream>

#include "Util.h"
#include "Args.h"
#include "Lineshape.h"
#include "JSON.h"

int main(int argc, char **argv) {
    Eigen::initParallel();

    args::ArgumentParser parser("Calculates different lineshapes then saves them as a splineshape for fast qMT lookup.\n"
                                "http://github.com/spinicist/QUIT");

    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag                    verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<std::string>  shape_arg(parser, "LINESHAPE", "Choose lineshape (Gauss/Lorentzian/Super-lorentzian)", {'l', "lineshape"}, "G");
    args::ValueFlag<double>       T2b(parser, "T2 BOUND", "Choose bound-pool T2 (default 10µs)", {'t', "T2b"}, 10e-6);
    args::ValueFlag<double>       frq_count(parser, "FREQUENCY COUNT", "Number of frequencies (default 10)", {'n', "frq_count"}, 10);
    args::ValueFlag<double>       frq_start(parser, "FREQUENCY START", "First saturation frequency (default 1000 Hz)", {'s', "frq_start"}, 1e3);
    args::ValueFlag<double>       frq_spacing(parser, "FREQUENCY SPACING", "Spacing of frequencies (default 1000 Hz)", {'p', "frq_space"}, 1e3);

    QI::ParseArgs(parser, argc, argv, verbose);
    QI_LOG(verbose, "Bound-pool T2:" << T2b.Get());
    auto frqs = Eigen::ArrayXd::LinSpaced(frq_count.Get(), frq_start.Get(), frq_start.Get() + frq_spacing.Get()*(frq_count.Get() - 1));
    QI_LOG(verbose, "Frequency count: " << frqs.rows());
    Eigen::ArrayXd values;
    if (shape_arg.Get() == "G") {
        QI::Lineshapes::Gaussian ls;
        values = ls.value(frqs, T2b.Get());
    } else {
        QI_FAIL("Unknown lineshape: " << shape_arg.Get());
    }
    const auto lineshape = QI::Lineshapes::Splineshape(frqs, values, T2b);
    rapidjson::Document doc;
    doc.SetObject();
    QI_JSONIFY(lineshape, doc);
    QI::WriteJSON(std::cout, doc);
    QI_LOG(verbose, "Finished.");
    return EXIT_SUCCESS;
}
