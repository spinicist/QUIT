#include "ModelHelpers.h"

namespace QI {

Eigen::ArrayXd add_noise(Eigen::ArrayXd const &s, double const sigma) {
    Eigen::ArrayXcd noise(s.rows());
    // Simple Box Muller transform
    Eigen::ArrayXd U = (Eigen::ArrayXd::Random(s.rows()) * 0.5) + 0.5;
    Eigen::ArrayXd V = (Eigen::ArrayXd::Random(s.rows()) * 0.5) + 0.5;
    noise.real()     = (sigma / M_SQRT2) * (-2. * U.log()).sqrt() * cos(2. * M_PI * V);
    noise.imag()     = (sigma / M_SQRT2) * (-2. * V.log()).sqrt() * sin(2. * M_PI * U);
    Eigen::ArrayXcd coutput(s.rows());
    coutput.real() = s + noise.real();
    coutput.imag() = noise.imag();
    return coutput.abs();
}

Eigen::ArrayXcd add_noise(Eigen::ArrayXcd const &s, double const sigma) {
    Eigen::ArrayXcd noise(s.rows());
    // Simple Box Muller transform
    Eigen::ArrayXd U = (Eigen::ArrayXd::Random(s.rows()) * 0.5) + 0.5;
    Eigen::ArrayXd V = (Eigen::ArrayXd::Random(s.rows()) * 0.5) + 0.5;
    noise.real()     = (sigma / M_SQRT2) * (-2. * U.log()).sqrt() * cos(2. * M_PI * V);
    noise.imag()     = (sigma / M_SQRT2) * (-2. * V.log()).sqrt() * sin(2. * M_PI * U);
    return s + noise;
}

} // namespace QI
