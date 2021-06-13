#include "OnePoolSignals.h"

namespace QI {

Eigen::ArrayXcd SPGR1(double const &          PD,
                      double const &          T1,
                      double const &          T2,
                      double const &          f0,
                      double const &          B1,
                      SPGREchoSequence const &s) {

    Eigen::ArrayXd const sa = sin(B1 * s.FA);
    Eigen::ArrayXd const ca = cos(B1 * s.FA);

    double     E1   = exp(-s.TR / T1);
    auto const echo = std::polar(PD * exp(-s.TE / T2), 2. * M_PI * f0 * s.TE);
    return echo * ((1. - E1) * sa) / (1. - E1 * ca);
}

Eigen::ArrayXcd SSFP1(double const &          PD,
                      double const &          T1,
                      double const &          T2,
                      double const &          f0,
                      double const &          B1,
                      const QI::SSFPSequence &s) {
    const double         E1     = exp(-s.TR / T1);
    const double         E2     = exp(-s.TR / T2);
    const Eigen::ArrayXd d      = (1 - E1 * cos(B1 * s.FA) - (E2 * E2) * (E1 - cos(B1 * s.FA)));
    const double         a      = E2;
    const Eigen::ArrayXd b      = E2 * (1 - E1) * (1 + cos(B1 * s.FA)) / d;
    const double         theta0 = 2.0 * M_PI * f0 * s.TR;
    const Eigen::ArrayXd theta  = theta0 + s.PhaseInc;

    Eigen::ArrayXcd Eth(theta.size());
    Eth.real()       = cos(theta);
    Eth.imag()       = sin(theta);
    const double psi = theta0 / 2.0;

    const Eigen::ArrayXcd G = std::polar(PD * sqrt(E2), psi) * (1 - E1) * sin(B1 * s.FA) / d;
    const Eigen::ArrayXcd M = G * (1. - a * Eth) / (1 - b * cos(theta));
    return M;
}

} // namespace QI
