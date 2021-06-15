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
    auto const echo = std::polar(exp(-s.TE / T2), 2. * M_PI * f0 * s.TE);
    return PD * echo * ((1. - E1) * sa) / (1. - E1 * ca);
}

Eigen::ArrayXcd SSFP1(double const &          PD,
                      double const &          T1,
                      double const &          T2,
                      double const &          f0,
                      double const &          B1,
                      const QI::SSFPSequence &s) {
    const double E1 = exp(-s.TR / T1);
    const double E2 = exp(-s.TR / T2);

    const double               psi   = 2. * M_PI * f0 * s.TR;
    const Eigen::ArrayXd       alpha = s.FA * B1;
    const Eigen::ArrayXd       theta = psi + s.PhaseInc;
    const Eigen::ArrayXd       d     = (1. - E1 * E2 * E2 - (E1 - E2 * E2) * cos(alpha));
    const std::complex<double> echo  = sqrt(E2) * std::polar(1., psi / 2.);
    const Eigen::ArrayXcd      G     = -PD * echo * (1 - E1) * sin(alpha) / d;
    const Eigen::ArrayXd       b     = E2 * (1. - E1) * (1. + cos(alpha)) / d;

    Eigen::ArrayXcd et(theta.size());
    et.real() = cos(-theta);
    et.imag() = sin(-theta);

    const Eigen::ArrayXcd M = G * (1. - E2 * et) / (1 - b * cos(theta));

    return M;
}

} // namespace QI
