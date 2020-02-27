/*
 *  rf_integral.cpp
 *
 *  Copyright (c) 2020 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

// #define QI_DEBUG_BUILD 1

#include "Args.h"
#include "JSON.h"
#include "Macro.h"
#include "Util.h"

#include "mupa_sequence.h"
#include "rufis_pulse.h"

double round_sig(double value, int digits) {
    if (value == 0.0) // otherwise it will return 'nan' due to the log10() of zero
        return 0.0;

    double factor = std::pow(10.0, digits - std::ceil(std::log10(std::abs(value))));
    return std::round(value * factor) / factor;
}

/*
 * Main
 */
int rf_sim_main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates the time integral of a pulse in JSON format "
                                "\nhttp://github.com/spinicist/QUIT");
    args::HelpFlag       help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more messages", {'v', "verbose"});
    args::Flag uT(parser, "uT", "Units are microTesla, not radians per second", {'u', "uT"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=hardware limit or $QUIT_THREADS)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());

    args::Positional<std::string> in_file(parser, "INPUT", "Input JSON file");
    // args::Positional<std::string> output_path(parser, "OUTPUT", "Output JSON file");

    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::CheckPos(in_file);

    double const R1 = 1. / 1.0;
    double const R2 = 1. / 0.1;
    double const PD = 1.;

    using AugMat = Eigen::Matrix<double, 4, 4>;
    using AugVec = Eigen::Vector<double, 4>;
    AugMat R;
    R << -R2, 0, 0, 0,      //
        0, -R2, 0, 0,       //
        0, 0, -R1, PD * R1, //
        0, 0, 0, 0;

    AugVec m0{0, 0, PD, 1.};

    QI::Log(verbose, "Reading pulses");
    json input = QI::ReadJSON(in_file.Get());

    std::vector<RFPulse>   input_pulses = input.at("pulses").get<std::vector<RFPulse>>();
    std::vector<PrepPulse> output_pulses;

    double const scale = uT ? 267.52219 : 1.0; // radians per second per uT

    for (auto const pulse : input_pulses) {
        double int_b1    = 0;
        double int_b1_sq = 0;
        double max_b1    = 0;
        double eff_tv    = 0;
        double eff_long  = 0;

        double tact_total = 0;
        AugMat C_rf       = AugMat::Identity();
        AugMat C_both     = AugMat::Identity();
        AugVec m_rf       = m0;
        for (long ii = 0; ii < pulse.B1x.rows(); ii++) {

            double const &B1x = pulse.B1x[ii] * scale;
            double const &B1y = pulse.B1y[ii] * scale;
            auto const    dt  = pulse.timestep[ii] * 1e-6;

            double b1_sq = (B1x * B1x + B1y * B1y);
            double b1    = sqrt(b1_sq);
            int_b1 += b1 * dt;
            int_b1_sq += b1_sq * dt;
            if (b1 > max_b1) {
                max_b1 = b1;
            }

            AugMat rf;
            rf << 0, 0, -B1y, 0, //
                0, 0, -B1x, 0,   //
                B1y, B1x, 0, 0,  //
                0, 0, 0, 0;

            AugMat const A_rf   = (rf * dt).exp();
            AugMat const both   = R + rf;
            AugMat const A_both = (both * dt).exp();

            C_rf   = A_rf * C_rf;
            C_both = A_both * C_both;
            m_rf   = A_rf * m_rf;
            if (b1_sq > 0.) {
                tact_total += dt;
            }
            eff_tv += sqrt(m_rf[0] * m_rf[0] + m_rf[1] * m_rf[1]) * dt;
            eff_long += std::abs(m_rf[2]) * dt;
        }

        double const eff_flip_onres = atan2(m_rf.head(2).norm(), m_rf[2]);
        // double const p1             = (int_b1 / tact_total) / max_b1;
        // double const p2             = (int_b1_sq / tact_total) / (max_b1 * max_b1);

        output_pulses.push_back(PrepPulse{round_sig(eff_flip_onres, 3),
                                          round_sig(int_b1_sq, 4),
                                          round_sig(eff_long, 4),
                                          round_sig(eff_tv, 4)});

        AugMat approx;
        approx << 0, 0, 0, 0, //
            0, 0, 0, 0,       //
            0, 0, exp(-R2 * eff_tv) * cos(eff_flip_onres),
            PD * (1 - exp(-R1 * eff_long)), //
            0, 0, 0, 1;
    }
    json output;
    output["pulses"] = output_pulses;
    fmt::print("{}\n", output.dump(2));
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
